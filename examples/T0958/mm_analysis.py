#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform

if(platform.system() == 'Darwin'):  # Mac system (local machine)
    OPENAWSEM_LOCATION = "/Users/mingchenchen/Documents/openmmawsem/openmmawsem/"
elif(platform.system() == 'Linux'):
    OPENAWSEM_LOCATION = '/projects/pw8/wl45/openmmawsem/'
else:
    print("system unknown")
sys.path.insert(0, OPENAWSEM_LOCATION)
from openmmawsem import *
from myFunctions_helper import *


parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False)
parser.add_argument("--crystal", action="store_true", default=False)
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--globular", action="store_true", default=False)
parser.add_argument("--hybrid", action="store_true", default=False)
parser.add_argument("-c", "--chain", type=str, default="A")

args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

proteinName = pdb_id = args.protein
chain=args.chain.upper()
pdb = f"{pdb_id}.pdb"

# print(args)
with open('analysis_commandline_args.txt', 'w') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

# for compute Q
input_pdb_filename, cleaned_pdb_filename = prepare_pdb("crystal_structure.pdb", chain)
ensure_atom_order(input_pdb_filename, quiet=1)

pdb_trajectory = read_trajectory_pdb_positions("movie.pdb")
#pdb_trajectory = read_trajectory_pdb_positions("openmmMovie.pdb")
oa = OpenMMAWSEMSystem("t0958-openmmawsem.pdb", chains=chain, k_awsem=1.0, xml_filename=OPENAWSEM_LOCATION+"awsem.xml") # k_awsem is an overall scaling factor that will affect the relevant temperature scales

# apply forces
# forceGroupTable_Rev = {11:"Con", 12:"Chain", 13:"Chi", 14:"Excluded", 15:"Rama", 16:"Direct",
#                   17:"Burial", 18:"Mediated", 19:"Fragment"}
forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
                    "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "beta_2":24,"beta_3":25,"pap":26, "Total":list(range(11, 32)),
                    "Water":[16, 18], "beta":[23, 24, 25], "pap":26, "Q":1}
#forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
#                    "Burial":17, "Mediated":18, "Contact":18, "Fragment":19, "Membrane":20, "ER":21,"TBM_Q":22, "beta_1":23, "Total":list(range(11, 26)),
#                    "Water":[16, 18], "beta":[23, 24, 25], "Q":1}

forces = [
    oa.q_value("crystal_structure-cleaned.pdb"),
    oa.con_term(),
    oa.chain_term(),
    oa.chi_term(),
    oa.excl_term(),
    oa.rama_term(),
    oa.rama_proline_term(),
    oa.rama_ssweight_term(),
    oa.contact_term(z_dependent=False),
    oa.apply_beta_term_1(),
    oa.apply_beta_term_2(),
    oa.apply_beta_term_3(),
    oa.pap_term(),
    oa.fragment_memory_term(frag_location_pre="./"),
    oa.er_term(),
    oa.tbm_q_term(k_tbm_q=2000),
    # oa.membrane_term(),
]
oa.addForcesWithDefaultForceGroup(forces)

# start simulation
collision_rate = 5.0 / picoseconds

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, Platform.getPlatformByName("OpenCL"))

#showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane","ER","TBM_Q","beta_1", "Total"]
showEnergy = ["Q", "Con", "Chain", "Chi", "Excluded", "Rama", "Contact", "Fragment", "Membrane","ER","TBM_Q","beta_1","beta_2","beta_3","pap", "Total"]
# print("Steps", *showEnergy)
print(" ".join(["{0:<8s}".format(i) for i in ["Steps"] + showEnergy]))
for step, pdb in enumerate(pdb_trajectory):
    simulation.context.setPositions(pdb.positions)
    e = []
    for term in showEnergy:
        if type(forceGroupTable[term]) == list:
            g = set(forceGroupTable[term])
        elif forceGroupTable[term] == -1:
            g = -1
        else:
            g = {forceGroupTable[term]}
        state = simulation.context.getState(getEnergy=True, groups=g)
        if term == "Q":
            termEnergy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        else:
            termEnergy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        e.append(termEnergy)
#     print(*e)
    print(" ".join([f"{step:<8}"] + ["{0:<8.2f}".format(i) for i in e]))
#         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
