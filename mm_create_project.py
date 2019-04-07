#!/usr/bin/env python3
import os
import argparse
import sys
import openmmawsem
import helperFunctions.myFunctions


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__author__ = 'Wei Lu'

parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("protein", help="The name of the protein, for example 1r69. \
            do: python3 ~/OPENAWSEM_LOCATION/mm_create_project.py")
parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False)
parser.add_argument("--extended", action="store_true", default=False)
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--hybrid", action="store_true", default=False)


args = parser.parse_args()

# Print if in debug
if args.debug:
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# Remove the pdb extension from the protein if it is there.
if len(args.protein) > 4 and args.protein[-4:] == '.pdb':
    name = args.protein[:-4]
else:
    name = args.protein
pdb = f"{name}.pdb"
chain = args.chain.upper()

# Log the command to a file
with open('create_project_commandline_args.txt', 'w') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

# If the file does not exist download it from the pdb
if not os.path.exists(f"crystal_structure.pdb"):
    pdb_list = [name]
    helperFunctions.myFunctions.downloadPdb(pdb_list)
    helperFunctions.myFunctions.cleanPdb(pdb_list, chain="-1", toFolder="cleaned_pdbs")
    do(f"cp cleaned_pdbs/{pdb} crystal_structure.pdb")

# If the chain is not specified then select all the chains
if chain == "-1":
    chain = helperFunctions.myFunctions.getAllChains("crystal_structure.pdb")
    print("Chains to simulate: ", chain)

# for compute Q
input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb("crystal_structure.pdb", chain)
openmmawsem.ensure_atom_order(input_pdb_filename)
# get fasta, pdb, seq file ready
openmmawsem.getSeqFromCleanPdb(input_pdb_filename, chains=chain, writeFastaFile=True)
do(f"cp crystal_structure.fasta {name}.fasta")

if args.extended:
    do(f"{__location__}/helperFunctions/fasta2pdb.py " + name)
    helperFunctions.myFunctions.add_chain_to_pymol_pdb(pdb)  # only work for one chain only now
else:
    do(f"cp crystal_structure.pdb {pdb}")

input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb(pdb, chain)
openmmawsem.ensure_atom_order(input_pdb_filename)

os.system(f"cp {__location__}/parameters/burial_gamma.dat .")
os.system(f"cp {__location__}/parameters/gamma.dat .")
os.system(f"cp {__location__}/parameters/membrane_gamma.dat .")
os.system(f"cp {__location__}/parameters/anti_* .")
os.system(f"cp {__location__}/parameters/para_* .")

do(f"python {__location__}/helperFunctions/Pdb2Gro.py crystal_structure.pdb amh-go.gro")

## ssweight
do("stride crystal_structure.pdb > ssweight.stride")
do(f"python {__location__}/helperFunctions/stride2ssweight.py > ssweight")

# below used for zim and zimPosition file
if args.membrane or args.hybrid:
    do("grep -E 'CB|CA  GLY' crystal_structure-cleaned.pdb > cbs.data")
    do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zimPosition""")
    do("python3 ~/opt/create_zim.py")


if args.frag:
    do(f"cp crystal_structure.fasta {name}.fasta")
    do(f"cp {__location__}/database/cullpdb_pc80_* .")
    do(f"python {__location__}/helperFunctions/MultCha_prepFrags_index.py \
    cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 1 9 > logfile" % name)
    helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")

do(f"cp {__location__}/mm_run.py .")
do(f"cp {__location__}/mm_analysis.py .")
do(f"cp {__location__}/params.py .")