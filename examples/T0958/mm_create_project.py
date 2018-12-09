#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
import platform

if(platform.system() == 'Darwin'):  # Mac system (local machine)
    OPENAWSEM_LOCATION = "/Users/mingchenchen/Documents/openmmawsem/openmmawsem"
elif(platform.system() == 'Linux'):
    OPENAWSEM_LOCATION = '/Users/mingchenchen/Documents/openmmawsem/openmmawsem'
else:
    print("system unknown")
sys.path.insert(0, OPENAWSEM_LOCATION)
from openmmawsem import *
from helperFunctions.myFunctions import *

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False)
parser.add_argument("--crystal", action="store_true", default=False)
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--hybrid", action="store_true", default=False)


args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

proteinName = pdb_id = args.protein
# chain='A'
# chain='ABC'
chain = args.chain.upper()

pdb = f"{pdb_id}.pdb"

# print(args)
with open('create_project_commandline_args.txt', 'w') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

# Download the file and rename it to crystal_structure.pdb
if not os.path.exists(f"crystal_structure.pdb"):
    downloadPdb([proteinName])
    do(f"cp original_pdbs/{pdb} crystal_structure.pdb")


if chain == "-1":
    chain = getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

# for compute Q
input_pdb_filename, cleaned_pdb_filename = prepare_pdb("crystal_structure.pdb", chain)
ensure_atom_order(input_pdb_filename)
# get fasta, pdb, seq file ready
getSeqFromCleanPdb(input_pdb_filename, chains=chain)

if not args.crystal:
    do(f"cp crystal_structure.fasta {pdb_id}.fasta")
    do("python ~/script/fasta2pdb.py "+proteinName)
    add_chain_to_pymol_pdb(pdb)  # only work for one chain only now
else:
    do(f"cp crystal_structure.pdb {pdb}")

input_pdb_filename, cleaned_pdb_filename = prepare_pdb(pdb, chain)
ensure_atom_order(input_pdb_filename)


#os.system(f"cp {OPENAWSEM_LOCATION}parameters/burial_gamma.dat .")
#os.system(f"cp {OPENAWSEM_LOCATION}parameters/gamma.dat .")
#os.system(f"cp {OPENAWSEM_LOCATION}parameters/membrane_gamma.dat .")

#do("python2 ~/opt/Pdb2Gro.py crystal_structure.pdb amh-go.gro")

## ssweight
#do("stride crystal_structure.pdb > ssweight.stride")
#do("python2 ~/opt/script/stride2ssweight.py > ssweight")

# below used for zim and zimPosition file
if args.membrane or args.hybrid:
    do("grep -E 'CB|CA  GLY' crystal_structure-cleaned.pdb > cbs.data")
    do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zimPosition""")
    do("python3 ~/opt/create_zim.py")


if args.frag:
    do(f"cp crystal_structure.fasta {proteinName}.fasta")
    do("cp ~/opt/database/cullpdb_pc80_* .")
    do("python2 ~/opt/script/MultCha_prepFrags_index.py \
    cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 1 9 > logfile" % proteinName)
    check_and_correct_fragment_memory("frags.mem")

do(f"cp {OPENAWSEM_LOCATION}mm_run.py .")
do(f"cp {OPENAWSEM_LOCATION}mm_analysis.py .")
