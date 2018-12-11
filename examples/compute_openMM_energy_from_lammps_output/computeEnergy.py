import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import numpy as np
import fileinput
from itertools import product
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
import seaborn as sns
from os import listdir

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata
import matplotlib as mpl

# from notebookFunctions import *
# from .. import notebookFunctions


import sys
try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.insert(0, OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *


simulation_platform = "OpenCL" # OpenCL, CUDA, CPU, or Reference
pdb_id = '1r69'
pdb = f"{pdb_id}.pdb"
chain='T'

# input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, chain)
input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"
# ensure_atom_order(input_pdb_filename)
# getSeqFromCleanPdb(input_pdb_filename, chains='A')


# ensure_atom_order(input_pdb_filename)
# a = open(pdb).read().split("END")
# os.system("rm openmmMovie.pdb")
# os.system(f"echo 'REMARK converted from awsem lammps output' >> openmmMovie.pdb")
# for i in range(30):
#     with open("tmp.pdb", "w") as out:
#         out.write(a[i])
#     input_pdb_filename, cleaned_pdb_filename =prepare_pdb("tmp.pdb", chain)
#     os.system(f"echo 'MODEL  {i+1}' >> openmmMovie.pdb")
#     os.system("cat tmp-openmmawsem.pdb >> openmmMovie.pdb")

pdb_trajectory = read_trajectory_pdb_positions("openmmMovie.pdb")
oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, xml_filename="../../awsem.xml") # k_awsem is an overall scaling factor that will affect the relevant temperature scales

# apply forces
# forceGroupTable_Rev = {11:"Con", 12:"Chain", 13:"Chi", 14:"Excluded", 15:"Rama", 16:"Direct",
#                   17:"Burial", 18:"Mediated", 19:"Fragment"}
forceGroupTable = {"Con":11, "Chain":12, "Chi":13, "Excluded":14, "Rama":15, "Direct":16,
                  "Burial":17, "Mediated":18, "Fragment":19, "Total":-1,
                  "Water":[16, 18]}
forces = [
    oa.con_term(),
    oa.chain_term(),
    oa.chi_term(),
    oa.excl_term(),
    oa.rama_term(),
    oa.rama_proline_term(),
    oa.direct_term(),
    oa.burial_term(),
    oa.mediated_term(),
    oa.fragment_memory_term(frag_location_pre="./")
]
oa.addForcesWithDefaultForceGroup(forces)

# start simulation
collision_rate = 5.0 / picoseconds

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, Platform.getPlatformByName("OpenCL"))

showEnergy = ["Con", "Chain", "Chi", "Excluded", "Rama", "Water", "Burial", "Fragment", "Total"]
# showEnergy = ["Con", "Chain", "Chi", "Excluded", "Rama", "Water", "Burial", "Total"]  # frag removed in this example.
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
        termEnergy = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        e.append(termEnergy)
#     print(*e)
    print(" ".join([f"{step:<8}"] + ["{0:<8.2f}".format(i) for i in e]))
#         print(forceGroupTable[term], state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))