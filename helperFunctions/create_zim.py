#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import myPersonalFunctions
import fileinput
from small_script.myFunctions import *

# parser = argparse.ArgumentParser(
#     description="The goal of this python3 code is to automatically create \
#     the project template as fast as possible. Written by Wei Lu."
# )
# parser.add_argument("protein", help="The name of the protein")
# parser.add_argument("-d", "--debug", action="store_true", default=False)
# parser.add_argument("--frag", action="store_true", default=False)
# parser.add_argument("--crystal", action="store_true", default=False)
# parser.add_argument("--membrane", action="store_true", default=False)
# parser.add_argument("--globular", action="store_true", default=False)
# parser.add_argument("--hybrid", action="store_true", default=False)


# args = parser.parse_args()

seq = ""
with open("crystal_structure.fasta", "r") as f:
    for line in f:
        if line[0] == ">":
            pass
        else:
            # print(line)
            seq += line.strip()
# create_zim(f"crystal.seq")
print(seq, len(seq))

data = read_hydrophobicity_scale(seq, isNew=False)
z = data["DGwoct"].values
np.savetxt("zim", z, fmt="%.2f")