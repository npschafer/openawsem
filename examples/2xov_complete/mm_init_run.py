#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
from time import sleep
import fileinput
import importlib.util

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulations")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("--name", default="simulation", help="Name of the simulation")
parser.add_argument("--to", default="./", help="location of movie file")
parser.add_argument("-c", "--chain", type=str, default="-1")
parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
parser.add_argument("-p", "--platform", type=str, default="OpenCL")
parser.add_argument("-s", "--steps", type=float, default=1e5, help="step size")
parser.add_argument("--tempStart", type=float, default=800, help="Starting temperature")
parser.add_argument("--tempEnd", type=float, default=200, help="Ending temperature")
parser.add_argument("--fromCheckPoint", type=str, default=None, help="The checkpoint file you want to start from")
parser.add_argument("-m", "--simulation_mode", type=int, default=1,
                help="default 1,\
                        0: constant temperature,\
                        1: temperature annealing")
parser.add_argument("--params", type=str, default="params.py")
parser.add_argument("--subMode", type=int, default=2)
parser.add_argument("-f", "--forces", default="forces_setup.py")
parser.add_argument("--parameters", default=None)
args = parser.parse_args()


do = os.system
cd = os.chdir

with open('commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')
print(' '.join(sys.argv))

# simulation_platform = "CPU"  # OpenCL, CUDA, CPU, or Reference
# simulation_platform = "OpenCL"
simulation_platform = args.platform
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    if args.thread != -1:
        platform.setPropertyDefaultValue("Threads", str(args.thread))
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")

# if mm_run.py is not at the same location of your setup folder.
setupFolderPath = os.path.dirname(args.protein)
proteinName = pdb_id = os.path.basename(args.protein)


pwd = os.getcwd()
toPath = os.path.abspath(args.to)
checkPointPath = None if args.fromCheckPoint is None else os.path.abspath(args.fromCheckPoint)
forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)
parametersLocation = "." if args.parameters is None else os.path.abspath(args.parameters)
os.chdir(setupFolderPath)



chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"

if chain == "-1":
    chain = getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

paramPath = args.params

if args.to != "./":
    # os.system(f"mkdir -p {args.to}")
    os.makedirs(toPath, exist_ok=True)
    os.system(f"cp {paramPath} {toPath}/params.py")
    os.system(f"cp crystal_structure.fasta {toPath}/")
    os.system(f"cp crystal_structure.pdb {toPath}/")
    # os.system(f"cp {pdb} {args.to}/{pdb}")
    # pdb = os.path.join(args.to, pdb)

# import args.params as params
spec = importlib.util.spec_from_file_location("params", paramPath)
params = importlib.util.module_from_spec(spec)
spec.loader.exec_module(params)

input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"


# start simulation
collision_rate = 5.0 / picoseconds
checkpoint_file = "checkpnt.chk"
checkpoint_reporter_frequency = 10000



snapShotCount = 400
stepsPerT = int(args.steps/snapShotCount)
Tstart = args.tempStart
Tend = args.tempEnd

if stepsPerT == 0:
    reporter_frequency = 4000
else:
    reporter_frequency = stepsPerT
# reporter_frequency = 4000

print(f"using force setup file from {forceSetupFile}")
spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
# print(spec)
forces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(forces)




oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=OPENAWSEM_LOCATION+"awsem.xml")  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
pull_part_forces = forces.pull_part(oa)
oa.addForces(pull_part_forces)
# integrator = LangevinIntegrator(800*kelvin, 1/picosecond, 5*femtoseconds)
integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, 10*femtoseconds)
# integrator = LangevinIntegrator(800*kelvin, 1/picosecond, 1*picoseconds)
# integrator = CustomIntegrator(0.001)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
simulation.context.setVelocitiesToTemperature(Tstart*kelvin)  # set the initial velocities of the atoms according to the desired starting temperature
# simulation.minimizeEnergy()  # first, minimize the energy to a local minimum to reduce any large forces that might be present

simulation.reporters.append(PDBReporter(os.path.join(toPath, "pulling.pdb"), 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))  # output energy and temperature during simulation
simulation.step(int(1e4))

simulation.saveState(os.path.join(toPath, 'output.xml'))