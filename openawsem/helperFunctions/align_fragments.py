#!/usr/bin/env python3
""" Align fragments from a fragment memory file to a fasta file. """
import argparse
from pathlib import Path
import pandas as pd
import shutil

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
    out=''
    len_previous_sequence=0
    for name, sequence in fasta_data.items():
        num_dashes_before = len_previous_sequence
        num_dashes_after = max_length - len(sequence) - num_dashes_before
        aligned_sequence = ('-' * num_dashes_before) + sequence + ('-' * num_dashes_after)
        out+=f">{name}\n{aligned_sequence}\n"
        len_previous_sequence+=len(sequence)
    return out

def create_staggered_fasta(fasta_path,frags_path, parent_dir, copy_path=None, output_path=None):
    fasta_path = Path(fasta_path)
    frags_path = Path(frags_path)
    parent_dir = Path(parent_dir)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Could not find fasta file: {fasta_path}")
    if not frags_path.exists():
        raise FileNotFoundError(f"Could not find fragment memory file: {frags_path}")
    if not parent_dir.exists():
        raise FileNotFoundError(f"Could not find: {parent_dir}.")
    if output_path is None:
        #Print to stdout
        output_path = Path('/dev/stdout')
    else:
        output_path = Path(output_path)
    # Parse the fasta file
    fasta_data = parse_fasta(fasta_path)
    
    with open(output_path, 'w') as aligned_fragments:
        # Write staggered sequences to output
        aligned_fragments.write(staggered_sequences_from_fasta(fasta_data))

        complete_sequence = ''.join(fasta_data.values())
        three_to_one = {
                        # Standard amino acids (common)
                        'ALA': 'A', 'LEU': 'L', 'GLY': 'G', 'VAL': 'V', 'ILE': 'I',
                        'SER': 'S', 'THR': 'T', 'CYS': 'C', 'PRO': 'P', 'ASP': 'D',
                        'GLU': 'E', 'PHE': 'F', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H',
                        'ASN': 'N', 'GLN': 'Q', 'MET': 'M', 'TRP': 'W', 'TYR': 'Y',
                        # Modified amino acids and non-standard amino acids (less common but notable)
                        'UNK': 'X',   # Unknown amino acid
                        'MSE': 'M',  # Selenomethionine, used frequently in X-ray crystallography
                        # 'SEC': 'C',  # (U) Selenocysteine, known as the 21st amino acid
                        # 'PYL': 'K',  # (O) Pyrrolysine, the 22nd amino acid
                        # 'PTR': 'Y',  # O-phosphotyrosine, important in signaling
                        # 'TPO': 'T',  # Phosphothreonine
                        # 'SEP': 'S',  # Phosphoserine
                        # 'PCA': 'Q',  # Pyroglutamic acid, a cyclized form of glutamine
                        # 'HYP': 'P',  # Hydroxyproline, common in collagen
                        'M3L': 'K',  # N3-methyllysine, common in histones
                        # 'ORN': 'K',  # Ornithine, used in some engineered proteins
                        # 'CIT': 'R',  # Citrulline, commonly studied in the urea cycle and immune responses
                        # 'HIC': 'R',  # Homoarginine, arises through post-translational modifications
                        # 'CYD': 'C',  # Dicysteine
                        'CAS': 'C',  # S-dicysteine
                        # 'NLE': 'L',  # Norleucine, used in protein engineering and as a methionine surrogate
                        # 'NVA': 'V',  # Norvaline, used in metabolic studies and protein structure analysis
                        # 'BMT': 'T',  # Beta-methylthreonine, found in some antibiotics
                        # 'DHA': 'S',  # Dehydroalanine, formed by post-translational modification
                        # '5HP': 'P',  # 5-hydroxyproline, another modification found in collagen
                        # 'SOC': 'C',  # Cysteic acid, an oxidation product of cysteine
                        # 'OCS': 'C',  # Cysteine sulfonic acid, another oxidation product
                        # 'MHO': 'M',  # Methionine sulfoxide, results from oxidative stress
                        # 'PLQ': 'W'   # Pyrroloquinoline quinone adduct, found in some enzyme studies
                    }
        
        fragment_table = pd.read_csv(frags_path, skiprows=4, delim_whitespace=True, names=['Gro', 'seq_init', 'mem_init', 'length', 'strength'])
        #strange_aminoacids = []
        for mem in fragment_table.itertuples():
            gro_file = Path(parent_dir) / mem.Gro
            if copy_path is not None:
                shutil.copy(gro_file, Path(copy_path))
            gro = pd.read_csv(gro_file, skiprows=2, delim_whitespace=True, names=['ResID', 'resName', 'name', 'serial', 'x', 'y', 'z'])
            #If not in the three_to_one dictionary, replace with 'X'
            #strange_aminoacids+= [aa for aa in gro['resName'].unique() if aa not in three_to_one.keys()]
            fragment = ''.join(gro[(gro['name'] == 'CA') & (gro['ResID'] >= mem.mem_init) & (gro['ResID'] < (mem.mem_init + mem.length))]['resName'].map(three_to_one).fillna('X'))
            if len(fragment) != mem.length:
                print(f"Fragment length mismatch for fragment {fragment} in {mem.Gro}: {len(fragment)} != {mem.length}\n{mem}")
            sequence = '-'*(mem.seq_init-1) + fragment + '-'*(len(complete_sequence) - mem.length - mem.seq_init)
            aligned_fragments.write(f">{Path(mem.Gro).stem}\n{sequence}\n")
        #if strange_aminoacids:
        #    print(f"Unrecognized amino acids: {set(strange_aminoacids)}")


def main(args=None):
    parser = argparse.ArgumentParser(description="Process fasta and fragment files for alignment.")
    parser.add_argument("--fasta", required=True, help="Path to the fasta file.")
    parser.add_argument("--frags", help="Path to the fragment memory file.",default='frags.mem')
    parser.add_argument("--working_directory", help="Parent directory for file processing.",default='.')
    parser.add_argument("--output", help="Output file path for aligned fragments.",default=None)
    parser.add_argument("--copy", help="Move fragments in the fragment memory fileto another directory",default=None)
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    create_staggered_fasta(fasta_path=args.fasta, frags_path=args.frags, parent_dir=args.working_directory, copy_path=args.copy, output_path=args.output),

if __name__ == "__main__":
    main()