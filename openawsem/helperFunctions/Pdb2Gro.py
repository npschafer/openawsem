#!python
# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian
#
# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/
#
# Last Update: 10/03/2023
# -------------------------------------------------------------------------

import argparse
from Bio.PDB.PDBParser import PDBParser


class Atom:
    atom_no = 0
    atom_name = ""
    res_no = 0
    res_name = ""
    x = 0.0
    y = 0.0
    z = 0.0

    def __init__(self, atom_no, atom_name, res_no, res_name, x, y, z, desc=''):
        self.atom_no = atom_no
        self.atom_name = atom_name
        self.res_no = res_no
        self.res_name = res_name
        self.x = x
        self.y = y
        self.z = z
        self.desc = desc

    def __init__(self, atom_no, atom_name, res_no, res_name, xyz, desc=''):
        self.atom_no = atom_no
        self.atom_name = atom_name
        self.res_no = res_no
        self.res_name = res_name
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.desc = desc

    def print_(self):
        print(self.atom_no, self.atom_name, self.res_no, self.res_name, self.x, self.y, self.z, self.desc)

    def write_(self, f):
        f.write(("     " + str(self.res_no))[-5:])
        f.write(("     " + self.res_name)[-5:])
        f.write(" " + (self.atom_name + "    ")[:4])
        f.write(("     " + str(self.atom_no))[-5:])
        f.write(("        " + str(round(self.x / 10, 3)))[-8:])
        f.write(("        " + str(round(self.y / 10, 3)))[-8:])
        f.write(("        " + str(round(self.z / 10, 3)))[-8:])
        f.write("\n")


def pdb2gro(pdb_file, output, chain_name=""):
    p = PDBParser(PERMISSIVE=1, QUIET=True)

    pdb_id = pdb_file[:-4] if pdb_file.lower().endswith(".pdb") else pdb_file
    s = p.get_structure(pdb_id, pdb_file)
    chains = s[0].get_list()

    atoms = []

    for chain in chains:
        if chain_name == "" or chain.get_id() == chain_name:
            ires = 0
            iatom = 0
            res_name = ""
            for res in chain:
                is_regular_res = res.has_id('N') and res.has_id('CA') and res.has_id('C')
                res_id = res.get_id()[0]
                if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
                    ires += 1
                    res_name = res.get_resname()
                    residue_no = res.get_id()[1]
                    for atom in res:
                        iatom += 1
                        atom_name = atom.get_name()
                        xyz = atom.get_coord()
                        atoms.append(Atom(iatom, atom_name, residue_no, res_name, xyz))

    with open(output, 'w') as out:
        out.write(" Structure-Based gro file\n")
        out.write(f"{len(atoms):12}\n")
        for iatom in atoms:
            iatom.write_(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PDB to GRO format.")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("output", help="Output GRO file")
    parser.add_argument("chain", nargs='?', default="", help="Optional chain name")
    args = parser.parse_args()

    pdb2gro(args.pdb_file, args.output, args.chain)