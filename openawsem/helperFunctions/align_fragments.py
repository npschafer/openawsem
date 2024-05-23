#!/usr/bin/env python3
""" Align fragments from a fragment memory file to a fasta file. """
import argparse
from pathlib import Path
import pandas as pd

def parse_fasta(filename):
    with open(filename, 'r') as file:
        fasta = {}
        sequence_name = None
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sequence_name = line[1:]
                if sequence_name in fasta:
                    raise ValueError(f"Duplicate sequence identifier found: {sequence_name}")
                fasta[sequence_name] = ""
            else:
                if sequence_name is None:
                    raise ValueError("File format error: sequence data without a header found.")
                fasta[sequence_name] += line.replace(" ", "")
        return fasta

def staggered_sequences_from_fasta(fasta_data):
    max_length = max(len(seq) for seq in fasta_data.values())
    out = ''
    for name, sequence in fasta_data.items():
        num_dashes_before = (max_length - len(sequence)) // 2
        num_dashes_after = max_length - len(sequence) - num_dashes_before
        aligned_sequence = ('-' * num_dashes_before) + sequence + ('-' * num_dashes_after)
        out += f">{name}\n{aligned_sequence}\n"
    return out

def create_staggered_fasta(fasta_path,frags_path,output_path):
    fasta_path = Path(fasta_path)
    frags_path = Path(frags_path)
    output_path = Path(output_path)

    # Parse the fasta file
    fasta_data = parse_fasta(fasta_path)
    
    with open(output_path, 'w+') as aligned_fragments:
        # Write staggered sequences to output
        aligned_fragments.write(staggered_sequences_from_fasta(fasta_data))

        complete_sequence = ''.join(fasta_data.values())
        three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        
        fragment_table = pd.read_csv(frags_path, skiprows=4, delim_whitespace=True, names=['Gro', 'seq_init', 'mem_init', 'length', 'strength'])
        for mem in fragment_table.itertuples():
            gro_file = Path(args.parent_dir) / mem.Gro
            gro = pd.read_csv(gro_file, skiprows=2, delim_whitespace=True, names=['ResID', 'resName', 'name', 'serial', 'x', 'y', 'z'])
            sequence = '-'*(mem.seq_init-1) + ''.join(gro[(gro['name'] == 'CA') & (gro['ResID'] >= mem.mem_init) & (gro['ResID'] <= (mem.mem_init + mem.length - 1))]['resName'].replace(three_to_one)) + '-'*(len(complete_sequence) - mem.length - mem.seq_init)
            aligned_fragments.write(f">{Path(mem.Gro).stem}\n{sequence}\n")


def main(args=None):
    parser = argparse.ArgumentParser(description="Process fasta and fragment files for alignment.")
    parser.add_argument("--fasta", required=True, help="Path to the fasta file.")
    parser.add_argument("--frags", required=True, help="Path to the fragment memory file.")
    parser.add_argument("--parent_dir", required=True, help="Parent directory for file processing.")
    parser.add_argument("--output", required=True, help="Output file path for aligned fragments.")
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    create_staggered_fasta(args.fasta, args.frags, args.output)

if __name__ == "__main__":
    main()