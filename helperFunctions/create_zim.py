#!python
import os
import argparse
import sys
from time import sleep
import subprocess
# import myPersonalFunctions
import fileinput
# from small_script.myFunctions import *

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

def read_hydrophobicity_scale(seq, isNew=False, tableLocation="~/openMM"):
    seq_dataFrame = pd.DataFrame({"oneLetterCode":list(seq)})
    # HFscales = pd.read_table("~/opt/small_script/Whole_residue_HFscales.txt")
    HFscales = pd.read_csv(f"{tableLocation}/helperFunctions/Whole_residue_HFscales.txt", sep="\t")
    if not isNew:
        # Octanol Scale
        # new and old difference is at HIS.
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS+" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    else:
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS0" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    HFscales_with_oneLetterCode = HFscales.assign(oneLetterCode=HFscales.AA.str.upper().map(code)).dropna()
    data = seq_dataFrame.merge(HFscales_with_oneLetterCode, on="oneLetterCode", how="left")
    return data


seq = ""
with open("crystal_structure.fasta", "r") as f:
    for line in f:
        if line[0] == ">":
            pass
        else:
            # print(line)
            seq += line.strip()
# create_zim(f"crystal.seq")
print((seq, len(seq)))

data = read_hydrophobicity_scale(seq, isNew=False)
z = data["DGwoct"].values
np.savetxt("zim", z, fmt="%.2f")