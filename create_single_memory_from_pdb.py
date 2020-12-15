#!/usr/bin/env python3
import os
import argparse
import sys
import openmmawsem
import helperFunctions.myFunctions

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__author__ = 'Wei Lu'

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create\
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("-p", "--protein", default="crystal_structure-cleaned.pdb", help="The name of the protein(1r69 for example)")
parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
args = parser.parse_args()


chain = args.chain
pdb = args.protein

if not os.path.exists(args.protein):
    print("ERROR: the pdb you specified is not exist")
    exit()
name = os.path.basename(args.protein)[:-4]

# If the chain is not specified then select all the chains
if chain == "-1":
    chain = helperFunctions.myFunctions.getAllChains(pdb, removeDNAchains=True)
    print("Chains info read from crystal_structure.pdb, chains to simulate: ", chain)


# get fasta, pdb, seq file ready
chain = helperFunctions.myFunctions.getAllChains(pdb)

for c in chain:
    # print(f"convert chain {c} of crystal structure to Gro file")
    os.system(f"python {__location__}/helperFunctions/Pdb2Gro.py {pdb} {name}_{c}.gro {c}")

seq_data = helperFunctions.myFunctions.seq_length_from_pdb(pdb, chain)
with open("single_frags.mem", "w") as out:
    out.write("[Target]\nquery\n\n[Memories]\n")
    for (chain_name, chain_start_residue_index, seq_length) in seq_data:
        # print(f"write chain {chain_name}")
        out.write(f"{name}_{chain_name}.gro {chain_start_residue_index} 1 {seq_length} 20\n")   # residue index in Gro always start at 1.