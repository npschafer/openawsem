#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform
import importlib.util

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

# __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
# __author__ = 'Wei Lu'

from openmmawsem import *
from helperFunctions.myFunctions import *


parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-c", "--chain", type=str, default="-1")
parser.add_argument("--thread", type=int, default=2, help="default is using 2 CPUs, -1 is using all")
parser.add_argument("-p", "--platform", type=str, default="CPU", help="Could be OpenCL, CUDA and CPU")
parser.add_argument("-t", "--trajectory", type=str, default="./movie.pdb")
parser.add_argument("--params", type=str, default=None)
parser.add_argument("-o", "--output", type=str, default=None, help="The location of file that show your energy and Q infomation.")
args = parser.parse_args()

if (args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

proteinName = pdb_id = args.protein
chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"

simulation_platform = args.platform
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    if args.thread != -1:
        platform.setPropertyDefaultValue("Threads", str(args.thread))
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")


# print(args)
with open('analysis_commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

if chain == "-1":
    chain = getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

# for compute Q
input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"

fileType = args.trajectory[-3:]
if fileType == "pdb":
    pdb_trajectory = md.load(args.trajectory, stride=1)
elif fileType == "dcd":
    pdb_trajectory = md.load(args.trajectory, top=input_pdb_filename, stride=1)
    # may use iterload if loading still too slow
else:
    print(f"Unknown fileType {fileType}")
# pdb_trajectory = read_trajectory_pdb_positions(args.trajectory)


# import args.params as params
# default is the params.py under the same folder of trajectory file
if args.params is None:
    location = os.path.dirname(args.trajectory)
    params_location = os.path.join(location, "params.py")
else:
    params_location = args.params
spec = importlib.util.spec_from_file_location("params", params_location)
params = importlib.util.module_from_spec(spec)
spec.loader.exec_module(params)


oa = OpenMMAWSEMSystem(input_pdb_filename, chains=chain, k_awsem=1.0, xml_filename=f"{OPENAWSEM_LOCATION}/awsem.xml")  # k_awsem is an overall scaling factor that will affect the relevant temperature scales

# apply forces
forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
                    "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "beta_2":24,"beta_3":25,"pap":26, "Total":list(range(11, 27)),
                    "Water":[16, 18], "Beta":[23, 24, 25], "Pap":26, "Rg_Bias":27, "Q":1, "Rg":2}
# forceGroupTable_Rev = {11:"Con", 12:"Chain", 13:"Chi", 14:"Excluded", 15:"Rama", 16:"Direct",
#                   17:"Burial", 18:"Mediated", 19:"Fragment"}
# forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
#                     "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "beta_2":24,"beta_3":25,"pap":26, "Total":list(range(11, 27)),
#                     "Water":[16, 18], "Beta":[23, 24, 25], "Pap":26, "Q":1}
# forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
#                    "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "Total":list(range(11, 26)),
#                    "Water":[16, 18], "beta":[23, 24, 25], "Q":1}
forces = [
    q_value(oa, "crystal_structure-cleaned.pdb"),
    con_term(oa),
    chain_term(oa),
    chi_term(oa),
    excl_term(oa, periodic=params.periodic),
    rama_term(oa),
    rama_proline_term(oa),
    rama_ssweight_term(oa),
    contact_term(oa, k_contact=params.k_contact, z_dependent=params.z_dependent, inMembrane=params.inMembrane,
                    k_relative_mem=params.k_relative_mem, periodic=params.periodic),
    beta_term_1(oa),
    beta_term_2(oa),
    beta_term_3(oa),
    pap_term_1(oa),
    pap_term_2(oa),
    fragment_memory_term(oa, frag_file_list_file="./frag.mem", UseSavedFragTable=True),
    # er_term(oa),
    # tbm_q_term(oa, k_tbm_q=2000),
    # membrane_term(oa, k_membrane=params.k_membrane, membrane_center=params.membrane_center),
    # rg_bias_term(oa, k_rg=params.k_rg, rg0=params.rg0),
    # rg_bias_term(oa, k_rg=params.k_rg, rg0=params.rg0, atomGroup=[1])
]

oa.addForcesWithDefaultForceGroup(forces)

# start simulation
collision_rate = 5.0 / picoseconds

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)

showEnergy = ["Q", "Rg", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane", "Beta", "Pap", "Rg_Bias", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane", "Beta", "Pap", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane","ER","TBM_Q","beta_1", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane","ER","TBM_Q","beta_1","beta_2","beta_3","pap", "Total"]
# print("Steps", *showEnergy)
if args.output is None:
    outFile = os.path.join(os.path.dirname(args.trajectory), "info.dat")
else:
    outFile = args.output
with open(outFile, "w") as out:
    line = " ".join(["{0:<8s}".format(i) for i in ["Steps"] + showEnergy])
    print(line)
    out.write(line+"\n")
    # for step, pdb in enumerate(pdb_trajectory):
    #     simulation.context.setPositions(pdb.positions)
    for step in range(len(pdb_trajectory)):
        simulation.context.setPositions(pdb_trajectory.openmm_positions(step))
        e = []
        for term in showEnergy:
            if type(forceGroupTable[term]) == list:
                g = set(forceGroupTable[term])
            elif forceGroupTable[term] == -1:
                g = -1
            else:
                g = {forceGroupTable[term]}
            state = simulation.context.getState(getEnergy=True, groups=g)
            if term == "Q" or term == "Rg":
                termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            else:
                termEnergy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            e.append(termEnergy)
    #     print(*e)
        line = " ".join([f"{step:<8}"] + ["{0:<8.2f}".format(i) for i in e])
        print(line)
        out.write(line+"\n")
    #         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
