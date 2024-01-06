from openmm.app import *
from openmm import *
from simtk.unit import * #simtk.unit
from openmm import app

#Loading Simulation Inputs
prmtop = app.AmberPrmtopFile("GPX6_wat.prmtop")
inpcrd = app.AmberInpcrdFile("GPX6_wat_init0.rst")

#Setting up the Simulation Environment

#OpenMM requires several pieces of information to get started, but assigning variables
#in stages can help keep things organized (and also easy to modify for future use.)

# Using Particle-mesh Ewald
system = prmtop.createSystem(nonbondedMethod = app.PME,
    # with a 10 Angstrom cutoff distance
    nonbondedCutoff = 1.0*nanometers,
    # keeping bonds between hydrogens and heavy atoms at a fixed length
    constraints = app.HBonds,
    # setting the error tolerance for the Ewald
    ewaldErrorTolerance = 0.001,
    # letting waters be flexible instead of solid triangles
    rigidWater = False,
    # maintains the center ochkf the periodic box at the center of mass
    removeCMMotion = True)

#you will also want to set up your thermostat to run in NVT (and add a barostat to run in NPT).
# Temperature, Friction Coefficient, Timestep
# This will be added to the larger simulation
thermostat = LangevinIntegrator(300*kelvin, 1/picoseconds, 0.001*picoseconds)

# Pressure, Temperature (only used for calculation),
# Frequency (how frequently the system should update the box size)
barostat = MonteCarloBarostat(1.0*atmosphere, 300.0*kelvin, 25)

# The barostat is added directly to the system rather than the larger
# simulation.
system.addForce(barostat)

#if you need distance or angle restraints between atoms, you can add them here.
#Add distance restraints
distance_restraints = HarmonicBondForce()
# First atom index (Python starts at zero!), Second atom index, Distance,
# Force Constant
# You can use as many addBond() function calls as you like, if there are
# multiple distances to restrain.
# But be careful not to duplicate any, as they won't overwrite each other,
# just additively stack up.
# distance_restraints.addBond(index1, index2, distance*angstroms,
#     force*kilocalories_per_mole/angstroms**2)
#
# # Add the collection of distance restraints to your system.
# system.addForce(distance_restraints)
#
# # Add angle restraints
# angle_restraints = HarmonicAngleForce()
# # First, second, and third atom indices, with the second atom as the vertex of
# # the angle, the desired angle in radians, and the force constant.
# # As above, add as many angle restraints as you need, but be careful not to
# # duplicate them.
# angle_restraints.addAngle(index1, index2, index3, angle,
#     force*kilocalories_per_mole/radians**2)
#
# # Add the collection of angle restraints to the system.
# system.addForce(angle_restraints)

#Now we need to set up the entire simulation with what we ve done so far.
# The topology of the system
sim = Simulation(prmtop.topology,
    # The system itself, with all the restraints, barostat, etc.
    system,
    # The thermostat, in this case the LangevinIntegrator we established before
    thermostat)

# Gets the starting position for every atom in the system and
# initializes the simulation

# Establishes the periodic boundary conditions of the system
# based on the box size.
sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()

# Filename to save the state data into.
sim.reporters.append(StateDataReporter("GPX6_output1.log",
    1000,   # number of steps between each save.
    step = True,             # writes the step number to each line
    time = True,             # writes the time (in ps)
    potentialEnergy = True,  # writes potential energy of the system (KJ/mole)
    kineticEnergy = True,    # writes the kinetic energy of the system (KJ/mole)
    totalEnergy = True,      # writes the total energy of the system (KJ/mole)
    temperature = True,      # writes the temperature (in K)
    volume = True,           # writes the volume (in nm^3)
    density = True))         # writes the density (in g/mL)

# Updates the checkpoint file every 10,000 steps
sim.reporters.append(CheckpointReporter('GPX6_chk_out.rst',20000))

# Adds another frame to the CHARMM-formatted DCD (which can be easily
# converted to .mdcrd by cpptraj) every 10,000 timesteps.
# Its good practice to keep the trajectory and checkpoint on the same
# write frequency, in case you need to stop a job and resume it later.
sim.reporters.append(DCDReporter("GPX6.dcd",20000))

# Running the Molecular Dynamics Simulation
# We have now set up the entire system for running molecular dynamics!
# But first, as is good practice, we want to minimize the system first to clear any potential clashes between atoms and residues.
# Minimize the system
sim.context.setPositions(inpcrd.getPositions(asNumpy=True))
sim.minimizeEnergy(maxIterations=2000)

#You can set the maxIterations to whatever value you want, but be careful not to minimize for too long,
# as you can wind up starting from a structure that may be very different from your original files.
#From here, you can begin the actual MD portion.

# This instructs the program to take 1,000,000 timesteps, which corresponds
# to 1.0 ns of MD based on a 2fs timestep.
sim.step(2000000)

#Restarting a Simulation from a Checkpoint

# Sometimes we may want to pick up a simulation where we left off.
# This requires only some slight modifications to the process.
# Most of the script you built above can be used, with a few exceptions.
# First, remove the following portion, as we no longer want the starting positions determined by the .inpcrd file.

# Gets the starting position for every atom in the system and
# initializes the simulation
#############################################
#sim.context.setPositions(inpcrd.getPositions(asNumpy=True))
#Next, add this section immediately before declaring the checkpoint file to use going forward.
#############################################
## Open checkpoint file
#with open("restart.chk",'rb') as f:
#    sim.context.loadCheckpoint(f.read())
##############################################
#crec que he de substituir la linia del inpcrd per aixo ultim.
