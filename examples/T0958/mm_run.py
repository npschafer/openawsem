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


if(platform.system() == 'Darwin'):  # Mac system (local machine)
    OPENAWSEM_LOCATION = "/Users/mingchenchen/Documents/openmmawsem/openmmawsem/"
elif(platform.system() == 'Linux'):
    OPENAWSEM_LOCATION = '/projects/pw8/wl45/openmmawsem/'
else:
    print("system unknown")
sys.path.insert(0, OPENAWSEM_LOCATION)
from openmmawsem import *
from helperFunctions.myFunctions import *

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulations")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("--name", default="simulation", help="Name of the simulation")
parser.add_argument("-c", "--chain", type=str, default="-1")
args = parser.parse_args()


do = os.system
cd = os.chdir

with open('commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


simulation_platform = "CPU"  # OpenCL, CUDA, CPU, or Reference
simulation_platform = "OpenCL"
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    platform.setPropertyDefaultValue("Threads", "1")
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")

proteinName = pdb_id = args.protein
chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"


if chain == "-1":
    chain = getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, chain)
ensure_atom_order(input_pdb_filename)
getSeqFromCleanPdb(input_pdb_filename, chains=chain)



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
    oa.er_term(),
    oa.tbm_q_term(k_tbm_q=2000),
    oa.apply_beta_term_1(),
    oa.apply_beta_term_2(),
    #oa.apply_beta_term_3(),
    oa.pap_term(),
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
simulation.reporters.append(PDBReporter("movie.pdb", reporter_frequency)) # output PDBs of simulated structures

print("Simulation Starts")
start_time = time.time()

for i in range(100):
    integrator.setTemperature(3*(200-i)*kelvin)
    simulation.step(10000)

#simulation.step(int(1e6))
simulation.reporters.append(CheckpointReporter(checkpoint_file, checkpoint_reporter_frequency)) # save progress during the simulation

time_taken = time.time() - start_time  # time_taken is in seconds
hours, rest = divmod(time_taken,3600)
minutes, seconds = divmod(rest, 60)
print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")
