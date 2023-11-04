#!/usr/bin/env python3
import os
import argparse
import openawsem.helperFunctions
from .Pdb2Gro import pdb2gro

__author__ = 'Wei Lu'

def create_single_memory(pdb, chain):
    if not os.path.exists(pdb):
        print("ERROR: the pdb you specified is not exist")
        exit()
    name = os.path.basename(pdb)[:-4]

    # If the chain is not specified then select all the chains
    if chain == "-1":
        chain = openawsem.helperFunctions.myFunctions.getAllChains(pdb, removeDNAchains=True)
        print("Chains info read from crystal_structure.pdb, chains to simulate: ", chain)


    # get fasta, pdb, seq file ready
    chain = openawsem.helperFunctions.myFunctions.getAllChains(pdb)

    for c in chain:
        # print(f"convert chain {c} of crystal structure to Gro file")
        
        #os.system(f"python {__location__}/helperFunctions/Pdb2Gro.py {pdb} {name}_{c}.gro {c}")
        pdb2gro(pdb,f'{name}_{c}.gro',c)

    seq_data = openawsem.helperFunctions.myFunctions.seq_length_from_pdb(pdb, chain)
    with open("single_frags.mem", "w") as out:
        out.write("[Target]\nquery\n\n[Memories]\n")
        for (chain_name, chain_start_residue_index, seq_length) in seq_data:
            # print(f"write chain {chain_name}")
            out.write(f"{name}_{chain_name}.gro {chain_start_residue_index} 1 {seq_length} 20\n")   # residue index in Gro always start at 1.
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description="The goal of this python3 code is to automatically create\
        the project template as fast as possible. Written by Wei Lu."
    )
    parser.add_argument("-p", "--protein", default="crystal_structure-cleaned.pdb", help="The name of the protein(1r69 for example)")
    parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
    args = parser.parse_args()
    chain = args.chain
    pdb = args.protein
    create_memory(pdb, chain)

