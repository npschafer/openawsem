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
parser.add_argument("protein", help="The name of the protein(1r69 for example): \
                        python3 ~/OPENAWSEM_LOCATION/mm_create_project.py 1r69 or you can specify the target pdb file you want\
                        to simulation. for example do: python3 ~/OPENAWSEM_LOCATION/mm_create_project.py YOUR_PDB_LOCATION/1r69.pdb")
parser.add_argument("-c", "--chain", default="-1", help="chains to be simulated, could be for example 'abc'.")
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("--frag", action="store_true", default=False, help="Generate fragment memories")
parser.add_argument("--extended", action="store_true", default=False, help="Start from extended structure")
parser.add_argument("--membrane", action="store_true", default=False)
parser.add_argument("--hybrid", action="store_true", default=False)
parser.add_argument("--verbose", action="store_true", default=False)

args = parser.parse_args()

# Print if in debug
if args.debug:
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# Log the command to a file
with open('create_project_commandline_args.txt', 'w') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

# if you provide the pdb then we will use it for the project(move to folder named original_pdbs to prevent from overwritten), otherwise we download it online.
if args.protein[-4:] == '.pdb':
    if not os.path.exists(args.protein):
        print("ERROR: the pdb you specified is not exist")
        exit()
    name = os.path.basename(args.protein)[:-4]
    pdb = os.path.basename(args.protein)
    do("mkdir -p original_pdbs")
    do(f"cp {args.protein} original_pdbs/")
elif args.protein[-6:] == ".fasta":
    print("Creating simulation folder from fasta file.")
    name = os.path.basename(args.protein)[:-6]
    # use pymol to generate the initial pdb.
    do(f"python3 {__location__}/helperFunctions/fasta2pdb.py {name} -f {args.protein}")
    helperFunctions.myFunctions.add_chain_to_pymol_pdb(f"{name}.pdb")  # only work for one chain only now
    do("mkdir -p original_fasta")
    do("mkdir -p original_pdbs")
    do(f"cp {args.protein} original_fasta/")
    do(f"cp {name}.pdb crystal_structure.pdb")  # treat as crystal_structure.
    pdb = f"{name}.pdb"
else:
    # If the file does not exist download it from the server
    name = args.protein
    pdb = f"{name}.pdb"
    pdb_list = [name]
    helperFunctions.myFunctions.downloadPdb(pdb_list)

if not os.path.exists(f"crystal_structure.pdb"):
    helperFunctions.myFunctions.cleanPdb([name], chain="-1", toFolder="cleaned_pdbs", verbose=args.verbose)
    do(f"cp cleaned_pdbs/{pdb} crystal_structure.pdb")

chain = args.chain.upper()
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
    do(f"python3 {__location__}/helperFunctions/fasta2pdb.py extended -f {name}.fasta")
    helperFunctions.myFunctions.add_chain_to_pymol_pdb("extended.pdb")  # only work for one chain only now
    input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb("extended.pdb", "A")
    openmmawsem.ensure_atom_order(input_pdb_filename)

do(f"cp crystal_structure.pdb {pdb}")
input_pdb_filename, cleaned_pdb_filename = openmmawsem.prepare_pdb(pdb, chain)
openmmawsem.ensure_atom_order(input_pdb_filename)

os.system(f"cp {__location__}/parameters/burial_gamma.dat .")
os.system(f"cp {__location__}/parameters/gamma.dat .")
os.system(f"cp {__location__}/parameters/membrane_gamma.dat .")
os.system(f"cp {__location__}/parameters/anti_* .")
os.system(f"cp {__location__}/parameters/para_* .")

do(f"python {__location__}/helperFunctions/Pdb2Gro.py crystal_structure.pdb {name}.gro")

## ssweight
do("stride crystal_structure.pdb > ssweight.stride")
do(f"python {__location__}/helperFunctions/stride2ssweight.py > ssweight")
protein_length = helperFunctions.myFunctions.getFromTerminal("wc ssweight").split()[0]
print(f"protein: {name}, length: {protein_length}")

seq_data = helperFunctions.myFunctions.seq_length_from_pdb("crystal_structure-cleaned.pdb", chain)
with open("single_frags.mem", "w") as out:
    out.write("[Target]\nquery\n\n[Memories]\n")
    for (chain_name, chain_start_residue_index, seq_length) in seq_data:
        out.write(f"{name}.gro {chain_start_residue_index} {chain_start_residue_index} {seq_length} 20\n")

# below used for zim and zimPosition file
if args.membrane or args.hybrid:
    do("grep -E 'CB|CA  GLY' crystal_structure-cleaned.pdb > cbs.data")
    do("""awk '{if($9>15) print "1"; else if($9<-15) print "3"; else print "2"}'  cbs.data  > zimPosition""")
    helperFunctions.myFunctions.create_zim(f"crystal_structure.fasta", tableLocation=f"{__location__}/helperFunctions")


if args.frag:
    do(f"cp crystal_structure.fasta {name}.fasta")
    do(f"cp {__location__}/database/cullpdb_pc80_* .")
    do(f"python {__location__}/helperFunctions/MultCha_prepFrags_index.py \
    cullpdb_pc80_res3.0_R1.0_d160504_chains29712 %s.fasta 20 1 9 > logfile" % name)
    helperFunctions.myFunctions.check_and_correct_fragment_memory("frags.mem")
    helperFunctions.myFunctions.relocate(fileLocation="frags.mem", toLocation="fraglib")
    # print(f"{__location__}//Gros/")
    helperFunctions.myFunctions.replace(f"frags.mem", f"{__location__}//Gros/", "./fraglib/")
    do("cp frags.mem frag_memory.mem")

do(f"cp {__location__}/mm_run.py .")
do(f"cp {__location__}/mm_analysis.py .")
# do(f"cp {__location__}/params.py .")
do(f"cp {__location__}/forces_setup.py .")

print(f"{args.protein} project folder created")