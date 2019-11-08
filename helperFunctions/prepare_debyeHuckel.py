#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform
import importlib.util
import numpy as np


parser = argparse.ArgumentParser(
    description="generate a array for debyeHuckel term."
)
parser.add_argument("fasta", help="The name of the fasta file")
parser.add_argument("-s", "--saveTo", default="charge.txt")
args = parser.parse_args()

# fastaFile = "/Users/weilu/Research/examples/openMM_simulation/1r69.fasta"
fastaFile = args.fasta
with open(fastaFile) as input_data:
    seq = ""
    for line in input_data:
        if(line[0] == ">"):
            # print(line)
            pass
        elif(line == "\n"):
            pass
        else:
            seq += line.strip("\n")
print(seq)

charge_list = []
for i, res in enumerate(seq):
    if res == "R" or res == "K":
        charge_list.append([i, 1.0])
    elif res == "D" or res == "E":
        charge_list.append([i, -1.0])
    else:
        charge_list.append([i, 0.0])

# saveTo = "/Users/weilu/Research/examples/openMM_simulation/1r69.charge"
saveTo = args.saveTo
print(f"write charge info to file: {saveTo}")
np.savetxt(saveTo, charge_list, fmt='%i %.1f')
