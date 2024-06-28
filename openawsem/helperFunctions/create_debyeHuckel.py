#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np

def generate_charge_array(fasta_file, save_to):
    with fasta_file.open() as input_data:
        seq = ""
        for line in input_data:
            if line.startswith(">"):
                pass
            elif line == "\n":
                pass
            else:
                seq += line.strip("\n")
    print(seq)

    charge_list = []
    for i, res in enumerate(seq):
        if res in ("R", "K"):
            charge_list.append([i, 1.0])
        elif res in ("D", "E"):
            charge_list.append([i, -1.0])
        else:
            charge_list.append([i, 0.0])

    print(f"write charge info to file: {save_to}")
    np.savetxt(save_to, charge_list, fmt='%i %.1f')

def main():
    parser = argparse.ArgumentParser(description="Generate an array for the Debye-Huckel term.")
    parser.add_argument("fasta", help="The name of the fasta file", type=Path)
    parser.add_argument("-s", "--saveTo", default="charge.txt", type=Path)
    args = parser.parse_args()

    generate_charge_array(args.fasta, args.saveTo)

if __name__ == "__main__":
    main()