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
parser.add_argument("-o", "--output", type=str, default=None, help="The Name of file that show your energy and Q infomation.")
parser.add_argument("--subMode", type=int, default=3)
parser.add_argument("-f", "--forces", default="forces_setup.py")
parser.add_argument("--parameters", default=None)
parser.add_argument("--fromOpenMMPDB", action="store_true", default=False)
parser.add_argument("--fasta", type=str, default="")

args = parser.parse_args()

with open('analysis_commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')
print(' '.join(sys.argv))

if (args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# if mm_run.py is not at the same location of your setup folder.
setupFolderPath = os.path.dirname(args.protein)
setupFolderPath = "." if setupFolderPath == "" else setupFolderPath
proteinName = pdb_id = os.path.basename(args.protein)
chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"

trajectoryPath = os.path.abspath(args.trajectory)
if args.output is None:
    outFile = os.path.join(os.path.dirname(trajectoryPath), "info.dat")
else:
    outFile = os.path.join(os.path.dirname(trajectoryPath), args.output)

forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)
parametersLocation = "." if args.parameters is None else os.path.abspath(args.parameters)
os.chdir(setupFolderPath)

simulation_platform = args.platform
platform = Platform.getPlatformByName(simulation_platform)
if simulation_platform == "CPU":
    if args.thread != -1:
        platform.setPropertyDefaultValue("Threads", str(args.thread))
    print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")



if chain == "-1":
    chain = getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

# for compute Q
if args.fromOpenMMPDB:
    input_pdb_filename = proteinName
    seq=read_fasta("crystal_structure.fasta")
    print(f"Using Seq:\n{seq}")
else:
    suffix = '-openmmawsem.pdb'
    if pdb_id[-len(suffix):] == suffix:
        input_pdb_filename = pdb_id
    else:
        input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"
    seq=None

if args.fasta == "":
    seq = None
else:
    seq = seq=read_fasta(args.fasta)
    print(f"Using Seq:\n{seq}")

fileType = trajectoryPath[-3:]
if fileType == "pdb":
    pdb_trajectory = md.load(trajectoryPath, stride=1)
elif fileType == "dcd":
    pdb_trajectory = md.load(trajectoryPath, top=input_pdb_filename, stride=1)
    # may use iterload if loading still too slow
else:
    print(f"Unknown fileType {fileType}")
# pdb_trajectory = read_trajectory_pdb_positions(trajectoryPath)



oa = OpenMMAWSEMSystem(input_pdb_filename, chains=chain, k_awsem=1.0, xml_filename=f"{OPENAWSEM_LOCATION}/awsem.xml", seqFromPdb=seq)  # k_awsem is an overall scaling factor that will affect the relevant temperature scales

print(f"using force setup file from {forceSetupFile}")
spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
# print(spec)
forces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(forces)
forces = forces.set_up_forces(oa, computeQ=True, submode=args.subMode, contactParameterLocation=parametersLocation)
oa.addForcesWithDefaultForceGroup(forces)
# print(forces)

# start simulation
collision_rate = 5.0 / picoseconds

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)

# apply forces
forceGroupTable = {"Backbone":20, "Rama":21, "Contact":22, "Fragment":23, "Membrane":24, "ER":25, "TBM_Q":26, "Beta":27, "Pap":28, "Helical":29,
                    "Q":1, "Rg":2, "Qc":3,
                    "Helix_orientation":18, "Pulling":19,
                    "Total":list(range(11, 32))
                    # , "Q_wat":4, "Q_mem":5, "Debye_huckel":30
                   }

# forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
#                     "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "beta_2":24,"beta_3":25,"pap":26, "Total":list(range(11, 32)),
#                     "Water":[16, 18], "Beta":[23, 24, 25], "Pap":26, "Rg_Bias":27, "Helical":28, "Pulling":29, "Q":1, "Rg":2, "Qc":3, "Q_wat":4, "Q_mem":5}
# forceGroupTable_Rev = {11:"Con", 12:"Chain", 13:"Chi", 14:"Excluded", 15:"Rama", 16:"Direct",
#                   17:"Burial", 18:"Mediated", 19:"Fragment"}
# forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
#                     "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "beta_2":24,"beta_3":25,"pap":26, "Total":list(range(11, 27)),
#                     "Water":[16, 18], "Beta":[23, 24, 25], "Pap":26, "Q":1}
# forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
#                    "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "Total":list(range(11, 26)),
#                    "Water":[16, 18], "beta":[23, 24, 25], "Q":1}

showEnergy = ["Q", "Rg", "Backbone", "Rama", "Contact", "Fragment", "Membrane", "ER", "TBM_Q", "Beta", "Pap", "Helical", "Total"]
# , "Disulfide"
# showEnergy = ["Q", "Qc", "Rg", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Helical", "Fragment", "Membrane", "ER", "Beta", "Pap", "Total"]
# showEnergy = ["Q", "Qc", "Rg", "Pulling", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "FamilyFold", "Fragment", "Membrane", "Beta", "Pap", "Rg_Bias", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane", "Beta", "Pap", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane","ER","TBM_Q","beta_1", "Total"]
# showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane","ER","TBM_Q","beta_1","beta_2","beta_3","pap", "Total"]
# print("Steps", *showEnergy)
print("Printing energies")

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
            if term == "Q" or term == "Rg" or term == "Qc" or term == "Q_wat" or term == "Q_mem":
                termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
            else:
                termEnergy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            e.append(termEnergy)
    #     print(*e)
        line = " ".join([f"{step:<8}"] + ["{0:<8.2f}".format(i) for i in e])
        print(line)
        out.write(line+"\n")
    #         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
