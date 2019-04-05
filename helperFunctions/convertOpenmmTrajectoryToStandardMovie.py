import os
import argparse
import sys
from myFunctions import *

parser = argparse.ArgumentParser(
    description="Convert openMM output to standard Pdbs."
)
parser.add_argument("openmm", help="The name of the OpenAWSEM output")
parser.add_argument("fasta", default="crystal_structure.fasta",help="Fasta file")
args = parser.parse_args()

movieFile = args.openmm

seq_dic = get_seq_dic(fasta=args.fasta)
convert_openMM_to_standard_pdb(fileName=movieFile, seq_dic=seq_dic)
