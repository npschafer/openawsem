import os
import argparse
import sys

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.insert(0, OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

# from myFunctions import *
from helperFunctions.myFunctions import *


parser = argparse.ArgumentParser(
    description="Convert openMM output to standard Pdbs."
)
parser.add_argument("openmm", help="The name of the OpenAWSEM output")
parser.add_argument("-f", "--fasta", default="../crystal_structure.fasta", help="Fasta file")
args = parser.parse_args()

movieFile = args.openmm

seq_dic = get_seq_dic(fasta=args.fasta)
convert_openMM_to_standard_pdb(fileName=movieFile, seq_dic=seq_dic, back=False)

os.system("cp ~/openmmawsem/helperFunctions/complete_2xov.tcl .")