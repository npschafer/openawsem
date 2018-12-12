#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
from time import sleep
import fileinput


try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.insert(0, OPENAWSEM_LOCATION)
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
parser.add_argument("--platform", type=str, default="OpenCL")
parser.add_argument("-s", "--steps", type=int, default=int(1e5), help="step size")
parser.add_argument("--simulation_mode", type=int, default=0,
                help="default 0: constant temperature,\
                        1: temperature annealing")
args = parser.parse_args()


do = os.system
cd = os.chdir

with open('commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


# simulation_platform = "CPU"  # OpenCL, CUDA, CPU, or Reference
# simulation_platform = "OpenCL"
simulation_platform = args.platform
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    if args.thread != -1:
        platform.setPropertyDefaultValue("Threads", str(args.thread))
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")

proteinName = pdb_id = args.protein
chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"


if chain == "-1":
    chain = getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

if args.to != "./":
    os.system(f"mkdir -p {args.to}")
    # os.system(f"cp {pdb} {args.to}/{pdb}")
    # pdb = os.path.join(args.to, pdb)

input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"

reporter_frequency = 10000
oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=OPENAWSEM_LOCATION+"awsem.xml") # k_awsem is an overall scaling factor that will affect the relevant temperature scales

# apply forces
forces = [
    oa.con_term(),
    oa.chain_term(),
    oa.chi_term(),
    oa.excl_term(),
    oa.rama_term(),
    oa.rama_proline_term(),
    oa.rama_ssweight_term(),
    oa.contact_term(z_dependent=False),
    # oa.er_term(),
    # oa.tbm_q_term(k_tbm_q=2000),
    # oa.apply_beta_term_1(),
    # oa.apply_beta_term_2(),
    # oa.apply_beta_term_3(),
    # oa.pap_term(),
    #oa.additive_amhgo_term(pdb_file = "1r69.pdb", chain_name="A"),
    #oa.direct_term(),
    #oa.burial_term(),
    #oa.mediated_term(),
    oa.fragment_memory_term(frag_location_pre="./"),
    #oa.membrane_term(),
]
oa.addForces(forces)

# start simulation
collision_rate = 5.0 / picoseconds
checkpoint_file = "restart"
checkpoint_reporter_frequency = 10000

integrator = LangevinIntegrator(600*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
simulation.context.setPositions(oa.pdb.positions) # set the initial positions of the atoms
# simulation.context.setVelocitiesToTemperature(300*kelvin) # set the initial velocities of the atoms according to the desired starting temperature
simulation.minimizeEnergy() # first, minimize the energy to a local minimum to reduce any large forces that might be present
simulation.reporters.append(StateDataReporter(stdout, reporter_frequency, step=True, potentialEnergy=True, temperature=True)) # output energy and temperature during simulation
simulation.reporters.append(PDBReporter(os.path.join(args.to, "movie.pdb"), reporter_frequency))  # output PDBs of simulated structures

print("Simulation Starts")
start_time = time.time()

if args.simulation_mode == 0:
    simulation.step(args.steps)
elif args.simulation_mode == 1:
    for i in range(100):
        integrator.setTemperature(3*(200-i)*kelvin)
        simulation.step(10000)

#simulation.step(int(1e6))
simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency)) # save progress during the simulation

time_taken = time.time() - start_time  # time_taken is in seconds
hours, rest = divmod(time_taken,3600)
minutes, seconds = divmod(rest, 60)
print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")
