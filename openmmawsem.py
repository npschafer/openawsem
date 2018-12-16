from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from pdbfixer import *
import mdtraj as md
from Bio.PDB.Polypeptide import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList
from Bio.PDB import PDBIO
from itertools import product, combinations
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import textwrap

se_map_3_letter = {'ALA': 0,  'PRO': 1,  'LYS': 2,  'ASN': 3,  'ARG': 4,
                   'PHE': 5,  'ASP': 6,  'GLN': 7,  'GLU': 8,  'GLY': 9,
                   'ILE': 10, 'HIS': 11, 'LEU': 12, 'CYS': 13, 'MET': 14,
                   'SER': 15, 'THR': 16, 'TYR': 17, 'VAL': 18, 'TRP': 19}

se_map_1_letter = {'A': 0,  'P': 1,  'K': 2,  'N': 3,  'R': 4,
                   'F': 5,  'D': 6,  'Q': 7,  'E': 8,  'G': 9,
                   'I': 10, 'H': 11, 'L': 12, 'C': 13, 'M': 14,
                   'S': 15, 'T': 16, 'Y': 17, 'V': 18, 'W': 19}


gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                            'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                            'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                            'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}

def identify_terminal_residues(pdb_filename):
    # identify terminal residues
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_filename)
    terminal_residues = {}
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            terminal_residues[chain.id] = (residues[0].id[1], residues[-1].id[1])
        return terminal_residues

def prepare_pdb(pdb_filename, chains_to_simulate):
    # for more information about PDB Fixer, see:
    # http://htmlpreview.github.io/?https://raw.github.com/pandegroup/pdbfixer/master/Manual.html
    # fix up input pdb
    cleaned_pdb_filename = "%s-cleaned.pdb" % pdb_filename[:-4]
    input_pdb_filename = "%s-openmmawsem.pdb" % pdb_filename[:-4]

    fixer = PDBFixer(filename=pdb_filename)

    # remove unwanted chains
    chains = list(fixer.topology.chains())
    chains_to_remove = [i for i, x in enumerate(chains) if x.id not in chains_to_simulate]
    fixer.removeChains(chains_to_remove)

    #Identify Missing Residues
    fixer.findMissingResidues()
    fixer.missingResidues = {}

    #Replace Nonstandard Residues
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    #Remove Heterogens
    fixer.removeHeterogens(False)

    #Add Missing Heavy Atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    #Add Missing Hydrogens
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(cleaned_pdb_filename, 'w'))

    #Read sequence
    structure = PDBParser().get_structure('X', cleaned_pdb_filename)

    # identify terminal residues
    terminal_residues = identify_terminal_residues(cleaned_pdb_filename)

    # process pdb for input into OpenMM
    #Selects only atoms needed for the awsem topology
    output = open(input_pdb_filename, 'w')
    counter=0
    for line in open(cleaned_pdb_filename):
        splitline = line.split()
        if len(line)>4 and line[0:4] == "ATOM":
            try:
                atom_index=line[6:11].strip()
                atom_type=line[12:16].strip()
                res_type=line[17:20].strip()
                chain=line[21].strip()
                res_index=line[22:26].strip()
                x=line[30:38].strip()
                y=line[38:46].strip()
                z=line[46:54].strip()
                element = line[76:78].strip()
            except ValueError:
                print(line)
                raise
        else:
            continue
        awsem_atoms = ["CA", "O", "CB", "C", "H", "N"]
        if int(res_index) == terminal_residues[chain][0]:
            awsem_atoms.remove("N")
            awsem_atoms.remove("H")
        if int(res_index) == terminal_residues[chain][1]:
            awsem_atoms.remove("C")

        if atom_type in awsem_atoms:
            line=list(line)
            if res_type == "GLY":
                line[17:20] = "IGL"
            elif res_type == "PRO":
                line[17:20] = "IPR"
            else:
                line[17:20] = "NGP"
            if atom_type == "CB":
                line[77] = "B"
            line=''.join(line)
            output.write(line)
            counter+=1
    #print("The system contains %i atoms"%counter)
    output.close()

    #Fix Virtual Site Coordinates:
    prepare_virtual_sites(input_pdb_filename)

    return input_pdb_filename, cleaned_pdb_filename

def prepare_virtual_sites(pdb_file):
    parser = PDBParser(QUIET=True)
    structure=parser.get_structure('X',pdb_file,)
    for model in structure:
        for chain in model:
            r_im={}
            r_i={}
            for residue in chain:
                r_im=r_i
                r_i={}
                for atom in residue:
                    r_i[atom.get_name()]=atom
                if 'N' in r_i:
                    r_i['N'].set_coord( 0.48318*r_im['CA'].get_coord()+ 0.70328*r_i['CA'].get_coord()- 0.18643 *r_im['O'].get_coord())
                if 'C' in r_im:
                    r_im['C'].set_coord(0.44365*r_im['CA'].get_coord()+ 0.23520*r_i['CA'].get_coord()+ 0.32115 *r_im['O'].get_coord())
                if 'H' in r_i:
                    r_i['H'].set_coord( 0.84100*r_im['CA'].get_coord()+ 0.89296*r_i['CA'].get_coord()- 0.73389 *r_im['O'].get_coord())
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

def build_lists_of_atoms(nres, residues):
    # build lists of atoms, residue types, and bonds
    atom_types=['n', 'h', 'ca', 'c', 'o', 'cb']
    res_types=[]
    atom_lists=dict(zip(atom_types,[[] for i in range(len(atom_types))]))
    for residue in residues:
        res_types.append(residue.name)
        atom_types=['n', 'h', 'ca', 'c', 'o', 'cb']
        residue_atoms = [x.index for x in residue.atoms()]
        if residue.index == 0:
            atom_types.remove('n')
            atom_lists['n'].append(-1)
            atom_types.remove('h')
            atom_lists['h'].append(-1)
        if residue.index == nres-1:
            atom_types.remove('c')
            atom_lists['c'].append(-1)
            pass
        if residue.name == "IGL" and 'cb' in atom_types:
            atom_types.remove('cb')
            atom_lists['cb'].append(-1)
        if residue.name in "IPR" and 'h' in atom_types:
            atom_types.remove('h')
            atom_lists['h'].append(-1)
        assert len(residue_atoms)==len(atom_types), '%s\n%s'%(str(residue_atoms),str(atom_types))
        atom_types=[a.name.lower() for a in residue._atoms] #Sometimes the atom order may be different
        for atom, atype in zip(residue_atoms, atom_types):
                atom_lists[atype].append(atom)

    return atom_lists, res_types

def build_lists_of_atoms_2(nres, residues, atoms):
    res_id = 0
    n = h = ca = c = o = cb = -1
    atom_types=['n', 'h', 'ca', 'c', 'o', 'cb']
    atom_types_table = {'N':'n', 'H':'h', 'CA':'ca', 'C':'c', 'O':'o', 'CB':'cb'}
    res_types=[]
    for residue in residues:
        res_types.append(residue.name)

    atom_lists=dict(zip(atom_types,[[-1]*nres for i in range(len(atom_types))]))
    for atom in atoms:
        atom_lists[atom_types_table[atom.name]][atom.residue.index] = atom.index

    return atom_lists, res_types


def ensure_atom_order(input_pdb_filename, quiet=1):
    # ensure order of ['n', 'h', 'ca', 'c', 'o', 'cb']
    # to be more specific, this ensure 'ca' always show up before 'c'.
    def first(t):
        return t[0]
    order_table = {'N':0, 'H':1, 'CA':2, 'C':3, 'O':4, 'CB':5}
    one_residue = []
    with open("tmp.pdb", "w") as out:
        with open(input_pdb_filename, "r") as f:
            all_lines = f.readlines()
            pre_res_id = all_lines[0].split()[5]
            for line in all_lines:
                info = line.split()
                if info[0]!="ATOM":
                    sorted_residue = sorted(one_residue, key=first)
                    for a in sorted_residue:
                        out.write(a[1])
                    one_residue = []
                    out.write(line)
                    continue
                res_id = info[5]
                if res_id != pre_res_id:
                    pre_res_id = res_id
                    # sort them in order
                    sorted_residue = sorted(one_residue, key=first)
                    for a in sorted_residue:
                        out.write(a[1])
                    if sorted_residue != one_residue:
                        if not quiet:
                            print("Reorder atom position")
                            print("Original, Changed to")
                        for i, t in zip(one_residue, sorted_residue):
                            if not quiet:
                                print(i[2], i[3], ", ", t[2], t[3])
                    one_residue = []
                atomType = info[2]
                one_residue.append((order_table[atomType], line, info[1], atomType))
    os.system(f"mv tmp.pdb {input_pdb_filename}")

def inWhichChain(residueId, chain_ends):
    chain_table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i, end_of_chain_resId in enumerate(chain_ends):
        if end_of_chain_resId < residueId:
            pass
        else:
            return chain_table[i]


def isChainEdge(residueId, chain_starts, chain_ends, n=2):
    # n is how far away from the two ends count as in chain edge.
    atBegin = False
    atEnd = False
    for i in range(n):
        if (residueId-i) in chain_starts:
            atBegin = True
    for i in range(n):
        if (residueId+i) in chain_ends:
            atEnd = True
    return (atBegin or atEnd)

def get_chain_starts_and_ends(all_res):
    chain_starts = []
    chain_ends = []
    pos = -1
    for i, res in enumerate(all_res):
        if res.chain.index != pos:
            chain_starts.append(i)
            pos = res.chain.index
            if i > 0:
                chain_ends.append(i-1)
    chain_ends.append(len(all_res)-1)
    return chain_starts, chain_ends

def setup_virtual_sites(nres, system, n, h, ca, c, o, cb, res_type, chain_starts, chain_ends):
    # set virtual sites
    for i in range(nres):
        if i not in chain_starts:
            n_virtual_site = ThreeParticleAverageSite(ca[i-1], ca[i], o[i-1],
                                                      0.48318, 0.70328, -0.18643)
            system.setVirtualSite(n[i], n_virtual_site)
            if not res_type[i] == "IPR":
                h_virtual_site = ThreeParticleAverageSite(ca[i-1], ca[i], o[i-1],
                                                          0.84100, 0.89296, -0.73389)
                system.setVirtualSite(h[i], h_virtual_site)
        if  i not in chain_ends:
            c_virtual_site = ThreeParticleAverageSite(ca[i], ca[i+1], o[i],
                                                      0.44365, 0.23520, 0.32115)
            # print("Virtual", c[i])
            system.setVirtualSite(c[i], c_virtual_site)

def setup_bonds(nres, n, h, ca, c, o, cb, res_type, chain_starts, chain_ends):
    bonds = []
    for i in range(nres):
        bonds.append((ca[i], o[i]))
        if not res_type[i] == "IGL":
            bonds.append((ca[i], cb[i]))
        if  i not in chain_ends:
            bonds.append((ca[i], ca[i+1]))
            bonds.append((o[i], ca[i+1]))

    for i in range(nres):
        if i not in chain_starts and not res_type[i] == "IGL":
            bonds.append((n[i], cb[i]))
        if i not in chain_ends and not res_type[i] == "IGL":
            bonds.append((c[i], cb[i]))
        if i not in chain_starts and i not in chain_ends:
            bonds.append((n[i], c[i]))
    return bonds

def read_fasta(fastaFile):
    with open(fastaFile) as input_data:
        data = ""
        for line in input_data:
            if(line[0] == ">"):
                print(line)
            elif(line == "\n"):
                pass
            else:
                data += line.strip("\n")
    return data

def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated


def getSeqFromCleanPdb(input_pdb_filename, chains='A', writeFastaFile=False):
    cleaned_pdb_filename = input_pdb_filename.replace("openmmawsem.pdb", "cleaned.pdb")
    pdb = input_pdb_filename.replace("-openmmawsem.pdb", "")
    fastaFile = pdb + ".fasta"
    ThreeToOne = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
        'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
        'TYR':'Y','VAL':'V'}

    s = PDBParser().get_structure("X", cleaned_pdb_filename)
    m = s[0] # model 0
    seq = ""
    if writeFastaFile:
        with open(fastaFile, "w") as out:
            for chain in chains:
                out.write(f">{pdb.upper()}:{chain.upper()}|PDBID|CHAIN|SEQUENCE\n")
                c = m[chain]
                chain_seq = ""
                for residue in c:
                    residue_name = residue.get_resname()
                    chain_seq += ThreeToOne[residue_name]
                out.write("\n".join(textwrap.wrap(chain_seq, width=80))+"\n")
                seq += chain_seq
    else:
        for chain in chains:
            c = m[chain]
            chain_seq = ""
            for residue in c:
                residue_name = residue.get_resname()
                chain_seq += ThreeToOne[residue_name]
            seq += chain_seq
    return seq

def fixPymolPdb(location):
    # location = "1r69.pdb"
    with open("tmp", "w") as out:
        with open(location, "r") as f:
            for line in f:
                info = list(line)
                if len(info) > 21:
                    info[21] = "A"
                out.write("".join(info))
    os.system(f"mv tmp {location}")

def download(pdb_id):
    if not os.path.isfile(f"{pdb_id}.pdb"):
        PDBList().retrieve_pdb_file(pdb_id.lower(), pdir='.', file_format='pdb')
        os.rename("pdb%s.ent" % pdb_id, f"{pdb_id}.pdb")

class OpenMMAWSEMSystem:
    def __init__(self, pdb_filename, chains='A', xml_filename='awsem.xml', k_awsem=1.0):
        # read PDB
        self.pdb = PDBFile(pdb_filename)
        self.forcefield = ForceField(xml_filename)
        self.system = self.forcefield.createSystem(self.pdb.topology)
        # define convenience variables
        self.nres = self.pdb.topology.getNumResidues()
        self.natoms = self.pdb.topology.getNumAtoms()
        self.residues = list(self.pdb.topology.residues())
        self.resi = [x.residue.index for x in list(self.pdb.topology.atoms())]
        # build lists of atoms and residue types
        # self.atom_lists,self.res_type=build_lists_of_atoms(self.nres, self.residues)
        self.atom_lists,self.res_type=build_lists_of_atoms_2(self.nres, self.residues, self.pdb.topology.atoms())
        # print(self.atom_lists,self.res_type)
        self.n =self.atom_lists['n']
        self.h =self.atom_lists['h']
        self.ca=self.atom_lists['ca']
        self.c =self.atom_lists['c']
        self.o =self.atom_lists['o']
        self.cb=self.atom_lists['cb']

        self.chain_starts, self.chain_ends = get_chain_starts_and_ends(self.residues)
        # setup virtual sites
        # if Segmentation fault in setup_virtual_sites, use ensure_atom_order
        setup_virtual_sites(self.nres, self.system, self.n, self.h, self.ca, self.c, self.o,
                            self.cb, self.res_type, self.chain_starts, self.chain_ends)
        # setup bonds
        self.bonds = setup_bonds(self.nres, self.n, self.h, self.ca, self.c, self.o,
                            self.cb, self.res_type, self.chain_starts, self.chain_ends)
        # identify terminal_residues
        # self.terminal_residues = identify_terminal_residues(pdb_filename)
        # set overall scaling
        self.k_awsem = k_awsem
        # keep track of force names for output purposes
        self.force_names = []
        # save seq info
        self.seq = getSeqFromCleanPdb(pdb_filename, chains=chains)


    def addForces(self, forces):
        for i, (force) in enumerate(forces):
            self.addForce(force)
            force.setForceGroup(i+1)

    def addForcesWithDefaultForceGroup(self, forces):
        for i, (force) in enumerate(forces):
            self.addForce(force)

    def read_reference_structure_for_q_calculation(self, pdb_file, chain_name, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.8*nanometers):
        structure_interactions = []
        parser = PDBParser()
        structure = parser.get_structure('X', pdb_file)
        chain = structure[0][chain_name]
        residues = [x for x in chain]
        for i, residue_i in enumerate(residues):
            for j, residue_j in enumerate(residues):
                ca_list = []
                cb_list = []
                atom_list_i = []
                atom_list_j = []
                if i-j >= min_seq_sep and i-j <= max_seq_sep:  # taking the signed value to avoid double counting
                    ca_i = residue_i['CA']
                    ca_list.append(ca_i)
                    atom_list_i.append(ca_i)
                    ca_j = residue_j['CA']
                    ca_list.append(ca_j)
                    atom_list_j.append(ca_j)
                    if not residue_i.get_resname() == "GLY":
                        cb_i = residue_i['CB']
                        cb_list.append(cb_i)
                        atom_list_i.append(cb_i)
                    if not residue_j.get_resname() == "GLY":
                        cb_j = residue_j['CB']
                        cb_list.append(cb_j)
                        atom_list_j.append(cb_j)
                    for atom_i, atom_j in product(atom_list_i, atom_list_j):
                        r_ijN = abs(atom_i - atom_j)/10.0*nanometers # convert to nm
                        if r_ijN <= contact_threshold:
                            sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                            gamma_ij = 1.0
                            if atom_i in ca_list:
                                i_index = self.ca[i]
                            if atom_i in cb_list:
                                i_index = self.cb[i]
                            if atom_j in ca_list:
                                j_index = self.ca[j]
                            if atom_j in cb_list:
                                j_index = self.cb[j]
                            structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                            structure_interactions.append(structure_interaction)

        return structure_interactions

    def read_reference_structure_for_q_calculation_2(self, pdb_file, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.8*nanometers):
        # default use all chains in pdb file.
        structure_interactions = []
        parser = PDBParser()
        structure = parser.get_structure('X', pdb_file)
        model = structure[0]
        chain_start = 0
        count = 0
        for chain in model.get_chains():
            chain_start += count
            count = 0
            for i, residue_i in enumerate(chain.get_residues()):
                count += 1
                #  print(i, residue_i)
                for j, residue_j in enumerate(chain.get_residues()):
                    ca_list = []
                    cb_list = []
                    atom_list_i = []
                    atom_list_j = []
                    if i-j >= min_seq_sep and i-j <= max_seq_sep:  # taking the signed value to avoid double counting
                        ca_i = residue_i['CA']
                        ca_list.append(ca_i)
                        atom_list_i.append(ca_i)
                        ca_j = residue_j['CA']
                        ca_list.append(ca_j)
                        atom_list_j.append(ca_j)
                        if not residue_i.get_resname() == "GLY":
                            cb_i = residue_i['CB']
                            cb_list.append(cb_i)
                            atom_list_i.append(cb_i)
                        if not residue_j.get_resname() == "GLY":
                            cb_j = residue_j['CB']
                            cb_list.append(cb_j)
                            atom_list_j.append(cb_j)
                        for atom_i, atom_j in product(atom_list_i, atom_list_j):
                            r_ijN = abs(atom_i - atom_j)/10.0*nanometers # convert to nm
                            if r_ijN <= contact_threshold:
                                sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                                gamma_ij = 1.0
                                if atom_i in ca_list:
                                    i_index = self.ca[i+chain_start]
                                if atom_i in cb_list:
                                    i_index = self.cb[i+chain_start]
                                if atom_j in ca_list:
                                    j_index = self.ca[j+chain_start]
                                if atom_j in cb_list:
                                    j_index = self.cb[j+chain_start]
                                structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                                structure_interactions.append(structure_interaction)

        return structure_interactions

    def read_reference_structure_for_q_calculation_3(self, pdb_file, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, Qflag = 0):
        # default use all chains in pdb file.
        # this change use the canonical Qw/Qo calculation for reference Q
        # for Qw calculation is 0; Qo is 1;
        structure_interactions = []
        parser = PDBParser()
        structure = parser.get_structure('X', pdb_file)
        model = structure[0]
        chain_start = 0
        count = 0;
        for chain in model.get_chains():
            chain_start += count
            count = 0
            for i, residue_i in enumerate(chain.get_residues()):
                #  print(i, residue_i)
                count +=1
                for j, residue_j in enumerate(chain.get_residues()):
                        if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                            ca_i = residue_i['CA']

                            ca_j = residue_j['CA']

                            r_ijN = abs(ca_i - ca_j)/10.0*nanometers # convert to nm
                            if Qflag ==1 and r_ijN >= contact_threshold: continue
                            sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                            gamma_ij = 1.0
                            i_index = self.ca[i+chain_start]
                            j_index = self.ca[j+chain_start]
                            structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                            structure_interactions.append(structure_interaction)
        return structure_interactions

    def read_reference_structure_for_q_calculation_4(self, contact_threshold,rnative_dat,  min_seq_sep=3, max_seq_sep=np.inf):
        # use contact matrix for Q calculation
        # this change use the canonical Qw/Qo calculation for reference Q
        # for Qw calculation is 0; Qo is 1;
        in_rnative = np.loadtxt(rnative_dat); ## read in rnative_dat file for Q calculation
        structure_interactions = []
        chain_start = 0
        count = 0;
        for i in range(self.nres):
            chain_start += count
            count = 0
            for j in range(self.nres):
                count +=1
                if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                    r_ijN = in_rnative[i][j]/10.0*nanometers # convert to nm
                    if r_ijN < contact_threshold: continue
                    sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    i_index = self.ca[i]
                    j_index = self.ca[j]
                    structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                    structure_interactions.append(structure_interaction)
        return structure_interactions

    def tbm_q_term(self, k_tbm_q, tbm_q_min_seq_sep = 2, tbm_q_cutoff=0.2*nanometers, tbm_q_well_width=0.1, target_q = 1.0):
        ### Added by Mingchen Chen
        ### this function is solely used for template based modelling from rnative.dat file
        ### for details, refer to Chen, Lin & Lu Wolynes JCTC 2018
        print("TBM_Q term ON");
        tbm_q = CustomCVForce("0.5*k_tbm_q*(q-q0)^2")
        q = self.q_value_dat(contact_threshold=tbm_q_cutoff,rnative_dat = "rnative.dat",  min_seq_sep=tbm_q_min_seq_sep, max_seq_sep=np.inf)
        tbm_q.addCollectiveVariable("q", q)
        tbm_q.addGlobalParameter("k_tbm_q", k_tbm_q)
        tbm_q.addGlobalParameter("q0", target_q)
        tbm_q.setForceGroup(22)
        return tbm_q

    def q_value(self, reference_pdb_file, reference_chain_name='A', min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers):
        ### Modified by Mingchen to compute canonical QW/QO
        # create bond force for q calculation
        qvalue = CustomBondForce("(1/normalization)*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
        qvalue.addPerBondParameter("gamma_ij")
        qvalue.addPerBondParameter("r_ijN")
        qvalue.addPerBondParameter("sigma_ij")
        # create bonds
        # structure_interactions = self.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
        structure_interactions = self.read_reference_structure_for_q_calculation_3(reference_pdb_file,
            min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=0)
        #print(len(structure_interactions))
        #print(structure_interactions)
        qvalue.addGlobalParameter("normalization", len(structure_interactions))
        for structure_interaction in structure_interactions:
            qvalue.addBond(*structure_interaction)
        qvalue.setForceGroup(1)
        return qvalue


    def q_value_dat(self,contact_threshold ,rnative_dat="rnative.dat",  min_seq_sep=3, max_seq_sep=np.inf):
        ### Added by Mingchen
        ### this function is solely used for template based modelling from rnative.dat file
        ### for details, refer to Chen, Lin & Lu Wolynes JCTC 2018
        qvalue_dat = CustomBondForce("(1/normalization)*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
        qvalue_dat.addPerBondParameter("gamma_ij")
        qvalue_dat.addPerBondParameter("r_ijN")
        qvalue_dat.addPerBondParameter("sigma_ij")
        structure_interactions_tbm_q = self.read_reference_structure_for_q_calculation_4(contact_threshold=contact_threshold,rnative_dat="rnative.dat",  min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep)
        qvalue_dat.addGlobalParameter("normalization", len(structure_interactions_tbm_q))
        for structure_interaction_tbm_q in structure_interactions_tbm_q:
            qvalue_dat.addBond(*structure_interaction_tbm_q)
        return qvalue_dat



    def addForce(self, force):
        self.system.addForce(force)

    def con_term(self, k_con=50208, bond_lengths=[.3816, .240, .276, .153]):
        # add con forces
        # 50208 = 60 * 2 * 4.184 * 100. kJ/nm^2, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_con *= self.k_awsem
        con = HarmonicBondForce()
        for i in range(self.nres):
            con.addBond(self.ca[i], self.o[i], bond_lengths[1], k_con)
            if not self.res_type[i] == "IGL":
                con.addBond(self.ca[i], self.cb[i], bond_lengths[3], k_con)
            if i not in self.chain_ends:
                con.addBond(self.ca[i], self.ca[i+1], bond_lengths[0], k_con)
                con.addBond(self.o[i], self.ca[i+1], bond_lengths[2], k_con)
        con.setForceGroup(11)   # start with 11, so that first 10 leave for user defined force.
        return con

    def chain_term(self, k_chain=50208, bond_lengths=[0.2459108, 0.2519591, 0.2466597]):
        # add chain forces
        # 50208 = 60 * 2 * 4.184 * 100. kJ/nm^2, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_chain *= self.k_awsem
        chain = HarmonicBondForce()
        for i in range(self.nres):
            if i not in self.chain_starts and not self.res_type[i] == "IGL":
                chain.addBond(self.n[i], self.cb[i], bond_lengths[0], k_chain)
            if i not in self.chain_ends and not self.res_type[i] == "IGL":
                chain.addBond(self.c[i], self.cb[i], bond_lengths[1], k_chain)
            if i not in self.chain_starts and i not in self.chain_ends:
                chain.addBond(self.n[i], self.c[i],  bond_lengths[2], k_chain)
        chain.setForceGroup(12)
        return chain

    def chi_term(self, k_chi=251.04, chi0=-0.71):
        # add chi forces
        # The sign of the equilibrium value is opposite and magnitude differs slightly
        # 251.04 = 60 * 4.184 kJ, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_chi *= self.k_awsem
        chi = CustomCompoundBondForce(4, "k_chi*(chi*norm-chi0)^2;"
                                         "chi=crossproduct_x*r_cacb_x+crossproduct_y*r_cacb_y+crossproduct_z*r_cacb_z;"
                                         "crossproduct_x=(u_y*v_z-u_z*v_y);"
                                         "crossproduct_y=(u_z*v_x-u_x*v_z);"
                                         "crossproduct_z=(u_x*v_y-u_y*v_x);"
                                         "norm=1/((u_x*u_x+u_y*u_y+u_z*u_z)*(v_x*v_x+v_y*v_y+v_z*v_z)*(r_cacb_x*r_cacb_x+r_cacb_y*r_cacb_y+r_cacb_z*r_cacb_z))^0.5;"
                                         "r_cacb_x=x1-x4;"
                                         "r_cacb_y=y1-y4;"
                                         "r_cacb_z=z1-z4;"
                                         "u_x=x1-x2; u_y=y1-y2; u_z=z1-z2;"
                                         "v_x=x3-x1; v_y=y3-y1; v_z=z3-z1;")
        chi.addGlobalParameter("k_chi", k_chi)
        chi.addGlobalParameter("chi0", chi0)
        for i in range(self.nres):
            if i not in self.chain_starts and i not in self.chain_ends and not self.res_type[i] == "IGL":
                chi.addBond([self.ca[i], self.c[i], self.n[i], self.cb[i]])
        chi.setForceGroup(13)
        return chi

    def excl_term(self, k_excl=8368, r_excl=0.35):
        # add excluded volume
        # Still need to add element specific parameters
        # 8368 = 20 * 4.184 * 100 kJ/nm^2, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_excl *= self.k_awsem
        excl = CustomNonbondedForce("k_excl*step(r0-r)*(r-r0)^2")
        excl.addGlobalParameter("k_excl", k_excl)
        excl.addGlobalParameter("r0", r_excl)
        for i in range(self.natoms):
            excl.addParticle()
        # print(self.ca)
        # print(self.bonds)
        # print(self.cb)
        excl.addInteractionGroup(self.ca, self.ca)
        excl.addInteractionGroup([x for x in self.cb if x > 0], [x for x in self.cb if x > 0])
        excl.addInteractionGroup(self.ca, [x for x in self.cb if x > 0])
        excl.addInteractionGroup(self.o, self.o)

        excl.setCutoffDistance(r_excl)
        excl.setNonbondedMethod(excl.CutoffNonPeriodic)
        excl.createExclusionsFromBonds(self.bonds, 1)
        excl.setForceGroup(14)
        return excl

    def rama_term(self, k_rama=8.368, num_rama_wells=3, w=[1.3149, 1.32016, 1.0264], sigma=[15.398, 49.0521, 49.0954], omega_phi=[0.15, 0.25, 0.65], phi_i=[-1.74, -1.265, 1.041], omega_psi=[0.65, 0.45, 0.25], psi_i=[2.138, -0.318, 0.78]):
        # add Rama potential
        # 8.368 = 2 * 4.184 kJ/mol, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_rama *= self.k_awsem
        rama_function = ''.join(["w%d*exp(-sigma%d*(omega_phi%d*phi_term%d^2+omega_psi%d*psi_term%d^2))+" \
                                % (i, i, i, i, i, i) for i in range(num_rama_wells)])[:-1]
        rama_function = '-k_rama*(' + rama_function + ");"
        rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                                psi_term{i}=cos(psi_{i}-psi0{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                                for i in range(num_rama_wells)])
        rama_string = rama_function+rama_parameters
        rama = CustomCompoundBondForce(5, rama_string)
        for i in range(num_rama_wells):
            rama.addGlobalParameter(f"k_rama", k_rama)
            rama.addGlobalParameter(f"w{i}", w[i])
            rama.addGlobalParameter(f"sigma{i}", sigma[i])
            rama.addGlobalParameter(f"omega_phi{i}", omega_phi[i])
            rama.addGlobalParameter(f"omega_psi{i}", omega_psi[i])
            rama.addGlobalParameter(f"phi0{i}", phi_i[i])
            rama.addGlobalParameter(f"psi0{i}", psi_i[i])
        for i in range(self.nres):
            if i not in self.chain_starts and i not in self.chain_ends and not self.res_type[i] == "IGL" and not self.res_type[i] == "IPR":
                rama.addBond([self.c[i-1], self.n[i], self.ca[i], self.c[i], self.n[i+1]])
        rama.setForceGroup(15)
        return rama

    def rama_proline_term(self, k_rama_proline=8.368, num_rama_proline_wells=2, w=[2.17, 2.15], sigma=[105.52, 109.09], omega_phi=[1.0, 1.0], phi_i=[-1.153, -0.95], omega_psi=[0.15, 0.15], psi_i=[2.4, -0.218]):
        # add Rama potential for prolines
        # 8.368 = 2 * 4.184 kJ/mol, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_rama_proline *= self.k_awsem
        rama_function = ''.join(["w_P%d*exp(-sigma_P%d*(omega_phi_P%d*phi_term%d^2+omega_psi_P%d*psi_term%d^2))+" \
                                % (i, i, i, i, i, i) for i in range(num_rama_proline_wells)])[:-1]
        rama_function = '-k_rama_proline*(' + rama_function + ");"
        rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0_P{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                                psi_term{i}=cos(psi_{i}-psi0_P{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                                for i in range(num_rama_proline_wells)])
        rama_string = rama_function+rama_parameters
        rama = CustomCompoundBondForce(5, rama_string)
        for i in range(num_rama_proline_wells):
            rama.addGlobalParameter(f"k_rama_proline", k_rama_proline)
            rama.addGlobalParameter(f"w_P{i}", w[i])
            rama.addGlobalParameter(f"sigma_P{i}", sigma[i])
            rama.addGlobalParameter(f"omega_phi_P{i}", omega_phi[i])
            rama.addGlobalParameter(f"omega_psi_P{i}", omega_psi[i])
            rama.addGlobalParameter(f"phi0_P{i}", phi_i[i])
            rama.addGlobalParameter(f"psi0_P{i}", psi_i[i])
        for i in range(self.nres):
            if i not in self.chain_starts and i not in self.chain_ends and self.res_type[i] == "IPR":
                rama.addBond([self.c[i-1], self.n[i], self.ca[i], self.c[i], self.n[i+1]])
        rama.setForceGroup(15)
        return rama

    def rama_ssweight_term(self, k_rama_ssweight=8.368, num_rama_wells=2, w=[2.0, 2.0],
                        sigma=[419.0, 15.398], omega_phi=[1.0, 1.0], phi_i=[-0.995, -2.25],
                        omega_psi=[1.0, 1.0], psi_i=[-0.82, 2.16], location_pre="./"):
        # add RamaSS potential
        # 8.368 = 2 * 4.184 kJ/mol, converted from default value in LAMMPS AWSEM
        # multiply interaction strength by overall scaling
        k_rama_ssweight *= self.k_awsem
        rama_function = ''.join(["wSS%d*ssweight(%d,resId)*exp(-sigmaSS%d*(omega_phiSS%d*phi_term%d^2+omega_psiSS%d*psi_term%d^2))+" \
                                % (i, i, i, i, i, i, i) for i in range(num_rama_wells)])[:-1]
        rama_function = '-k_rama_ssweight*(' + rama_function + ");"
        rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0SS{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                                psi_term{i}=cos(psi_{i}-psi0SS{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                                for i in range(num_rama_wells)])
        rama_string = rama_function+rama_parameters
        ramaSS = CustomCompoundBondForce(5, rama_string)
        ramaSS.addPerBondParameter("resId")
        for i in range(num_rama_wells):
            ramaSS.addGlobalParameter(f"k_rama_ssweight", k_rama_ssweight)
            ramaSS.addGlobalParameter(f"wSS{i}", w[i])
            ramaSS.addGlobalParameter(f"sigmaSS{i}", sigma[i])
            ramaSS.addGlobalParameter(f"omega_phiSS{i}", omega_phi[i])
            ramaSS.addGlobalParameter(f"omega_psiSS{i}", omega_psi[i])
            ramaSS.addGlobalParameter(f"phi0SS{i}", phi_i[i])
            ramaSS.addGlobalParameter(f"psi0SS{i}", psi_i[i])
        for i in range(self.nres):
            if i not in self.chain_starts and i not in self.chain_ends and not self.res_type[i] == "IGL" and not self.res_type == "IPR":
                ramaSS.addBond([self.c[i-1], self.n[i], self.ca[i], self.c[i], self.n[i+1]], [i])
        ssweight = np.loadtxt(location_pre+"ssweight")
        ramaSS.addTabulatedFunction("ssweight", Discrete2DFunction(2, self.nres, ssweight.flatten()))
        ramaSS.setForceGroup(15)
        return ramaSS

    def direct_term(self, k_direct=4.184*1.5):
        k_direct *= self.k_awsem
        # print(self.ca, self.cb)
        # print(self.bonds)
        # print(self.nres)  # print 181 for 2xov
        # print(self.resi)  # print the rsidues index for each atom
        cb = self.cb
        # gamma = 1
        r_min = .45
        r_max = .65
        eta = 50  # eta actually has unit of nm^-1.
        min_sequence_separation = 10  # means j-i > 9
        nwell = 1
        gamma_ijm = np.zeros((nwell, 20, 20))
        # read in seq data.
        seq = self.seq
        # read in gamma info
        gamma_direct, gamma_mediated = read_gamma("gamma.dat")

        direct = CustomNonbondedForce(f"-k_direct*gamma_ijm(0, resName1, resName2)*theta; \
        theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r))); \
        eta={eta}")
        # direct = CustomNonbondedForce(f"-k_direct;")
        # direct = CustomNonbondedForce(f"-k_direct*gamma_ijm(0, resName1, resName2);")
        # direct = CustomNonbondedForce(f"-k_direct*gamma_ijm(0, resName1, resName2)*r;")
        direct.addGlobalParameter("k_direct", k_direct)
        direct.addGlobalParameter("rmin", r_min)
        direct.addGlobalParameter("rmax", r_max)



        # add per-particle parameters
        direct.addPerParticleParameter("resName")

        for i in range(self.natoms):
            direct.addParticle([gamma_se_map_1_letter[seq[self.resi[i]]]])


        for m in range(nwell):
            count = 0
            for i in range(20):
                for j in range(i, 20):
                    gamma_ijm[m][i][j] = gamma_direct[count][0]
                    gamma_ijm[m][j][i] = gamma_direct[count][0]
                    count += 1

        direct.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, 20, 20, gamma_ijm.flatten()))


        # direct.addInteractionGroup([x for x in cb if x > 0], [x for x in cb if x > 0])
        # direct.addInteractionGroup([x if x > 0 else y for x,y in zip(cb,self.ca)], [x if x > 0 else y for x,y in zip(cb,self.ca)])
        # direct.createExclusionsFromBonds(self.bonds, 11)
        # replace cb with ca for GLY
        cb_fixed = [x if x > 0 else y for x,y in zip(cb,self.ca)]
        # add interaction that are cutoff away
        # don't use this for multi chain simulation.
        for i, x in enumerate(cb_fixed):
            # print(i, x)
            direct.addInteractionGroup([x], cb_fixed[i+min_sequence_separation:])
        # print(cb)

        direct.setForceGroup(16)
        return direct


    def burial_term(self, k_burial=4.184, fastaFile="FastaFileMissing"):
        k_burial *= self.k_awsem
        burial_kappa = 4.0
        burial_ro_min = [0.0, 3.0, 6.0]
        burial_ro_max = [3.0, 6.0, 9.0]
        seq = self.seq
        eta = 50  # eta actually has unit of nm^-1.
        r_min = .45
        r_max = .65
        burial_gamma = np.loadtxt("burial_gamma.dat")

        # return burial
        # if ( lc->chain_no[i]!=lc->chain_no[j] || abs(lc->res_no[j] - lc->res_no[i])>1 )
        burial = CustomGBForce()

        burial_gamma_ij = np.zeros((20, 3))
        burial.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(20, 3, burial_gamma.T.flatten()))

        burial.addPerParticleParameter("resName")
        burial.addPerParticleParameter("resId")
        burial.addPerParticleParameter("isCb")
        burial.addGlobalParameter("k_burial", k_burial)
        burial.addGlobalParameter("eta", eta)
        burial.addGlobalParameter("burial_kappa", burial_kappa)
        burial.addGlobalParameter("rmin", r_min)
        burial.addGlobalParameter("rmax", r_max)
        index = burial.addComputedValue("rho", "step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)
        # print(burial.getComputedValueParameters(index))

        # replace cb with ca for GLY
        cb_fixed = [x if x > 0 else y for x,y in zip(self.cb,self.ca)]
        none_cb_fixed = [i for i in range(self.natoms) if i not in cb_fixed]
        for i in range(self.natoms):
            burial.addParticle([gamma_se_map_1_letter[seq[self.resi[i]]], self.resi[i], int(i in cb_fixed)])
        for i in range(3):
            burial.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            burial.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            burial.addEnergyTerm(f"-0.5*isCb*k_burial*burial_gamma_ij(resName, {i})*\
                                        (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                        tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

        # burial.addEnergyTerm("-k_burial*rho", CustomGBForce.SingleParticle)
        # burial.addEnergyTerm("-k_burial", CustomGBForce.SingleParticle)


        # print(len(none_cb_fixed), len(cb_fixed))
        for e1 in none_cb_fixed:
            for e2 in none_cb_fixed:
                if e1 > e2:
                    continue
                burial.addExclusion(e1, e2)
        for e1 in none_cb_fixed:
            for e2 in cb_fixed:
                burial.addExclusion(e1, e2)

        burial.setForceGroup(17)
        return burial


    def mediated_term(self, k_mediated=4.184*1.5):
        k_mediated *= self.k_awsem
        # print(self.nres)  # print 181 for 2xov
        # print(self.resi)  # print the rsidues index for each atom
        # gamma = 1
        r_min = .45
        r_max = .65
        r_minII = .65
        r_maxII = .95
        eta = 50  # eta actually has unit of nm^-1.
        eta_sigma = 7.0
        rho_0 = 2.6
        min_sequence_separation = 10  # means j-i > 9
        nwell = 1
        water_gamma_ijm = np.zeros((nwell, 20, 20))
        protein_gamma_ijm = np.zeros((nwell, 20, 20))
        # read in seq data.
        seq = self.seq
        # read in gamma info
        gamma_direct, gamma_mediated = read_gamma("gamma.dat")

        # mediated = CustomNonbondedForce(f"-k_mediated*densityGamma*theta2; \
        # densityGamma=sigmawater_gamma_ijm(0, resName1, resName2); \
        # theta2=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r))); \
        # eta={eta}")
        # mediated = CustomNonbondedForce(f"rho;")

        mediated = CustomGBForce()

        for m in range(nwell):
            count = 0
            for i in range(20):
                for j in range(i, 20):
                    water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
                    water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
                    count += 1

        for m in range(nwell):
            count = 0
            for i in range(20):
                for j in range(i, 20):
                    protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
                    protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
                    count += 1
        mediated.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, 20, 20, water_gamma_ijm.flatten()))
        mediated.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, 20, 20, protein_gamma_ijm.flatten()))

        # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
        res_table = np.zeros((self.nres, self.nres))
        for i in range(self.nres):
            for j in range(self.nres):
                resId1 = i
                chain1 = inWhichChain(resId1, self.chain_ends)
                resId2 = j
                chain2 = inWhichChain(resId2, self.chain_ends)
                if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                    res_table[i][j] = 1
                else:
                    res_table[i][j] = 0
        mediated.addTabulatedFunction("res_table", Discrete2DFunction(self.nres, self.nres, res_table.T.flatten()))
        mediated.addPerParticleParameter("resName")
        mediated.addPerParticleParameter("resId")
        mediated.addPerParticleParameter("isCb")
        mediated.addGlobalParameter("k_mediated", k_mediated)
        mediated.addGlobalParameter("eta", eta)
        mediated.addGlobalParameter("eta_sigma", eta_sigma)
        mediated.addGlobalParameter("rho_0", rho_0)
        mediated.addGlobalParameter("min_sequence_separation", min_sequence_separation)
        mediated.addGlobalParameter("rmin", r_min)
        mediated.addGlobalParameter("rmax", r_max)
        mediated.addGlobalParameter("rminII", r_minII)
        mediated.addGlobalParameter("rmaxII", r_maxII)

        mediated.addComputedValue("rho", "step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)
        # print(burial.getComputedValueParameters(index))

        # replace cb with ca for GLY
        cb_fixed = [x if x > 0 else y for x,y in zip(self.cb,self.ca)]
        none_cb_fixed = [i for i in range(self.natoms) if i not in cb_fixed]
        for i in range(self.natoms):
            mediated.addParticle([gamma_se_map_1_letter[seq[self.resi[i]]], self.resi[i], int(i in cb_fixed)])

        mediated.addEnergyTerm("-res_table(resId1, resId2)*k_mediated*thetaII*\
                                (sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(0, resName1, resName2));\
                                sigma_protein=1-sigma_water;\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)
        # print(len(none_cb_fixed), len(cb_fixed))
        for e1 in none_cb_fixed:
            for e2 in none_cb_fixed:
                if e1 > e2:
                    continue
                mediated.addExclusion(e1, e2)
        for e1 in none_cb_fixed:
            for e2 in cb_fixed:
                mediated.addExclusion(e1, e2)

        mediated.setForceGroup(18)
        return mediated


    def contact_term(self, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False):
        k_contact *= self.k_awsem
        # combine direct, burial, mediated.
        # default membrane thickness 1.5 nm

        r_min = .45
        r_max = .65
        r_minII = .65
        r_maxII = .95
        eta = 50  # eta actually has unit of nm^-1.
        eta_sigma = 7.0
        rho_0 = 2.6
        min_sequence_separation = 10  # means j-i > 9
        min_sequence_separation_mem = 13
        nwell = 2
        eta_switching = 10
        gamma_ijm = np.zeros((nwell, 20, 20))
        water_gamma_ijm = np.zeros((nwell, 20, 20))
        protein_gamma_ijm = np.zeros((nwell, 20, 20))

        # read in seq data.
        seq = self.seq
        # read in gamma info
        gamma_direct, gamma_mediated = read_gamma("gamma.dat")

        burial_kappa = 4.0
        burial_ro_min = [0.0, 3.0, 6.0]
        burial_ro_max = [3.0, 6.0, 9.0]
        burial_gamma = np.loadtxt("burial_gamma.dat")

        k_relative_mem = 400  # adjust the relative strength of gamma
        inMembrane = int(inMembrane)
        contact = CustomGBForce()

        m = 0  # water environment
        count = 0
        for i in range(20):
            for j in range(i, 20):
                gamma_ijm[m][i][j] = gamma_direct[count][0]
                gamma_ijm[m][j][i] = gamma_direct[count][0]
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
                water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
                protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
                count += 1
        # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
        res_table = np.zeros((nwell, self.nres, self.nres))
        for i in range(self.nres):
            for j in range(self.nres):
                resId1 = i
                chain1 = inWhichChain(resId1, self.chain_ends)
                resId2 = j
                chain2 = inWhichChain(resId2, self.chain_ends)
                if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                    res_table[0][i][j] = 1
                else:
                    res_table[0][i][j] = 0


        if z_dependent or inMembrane:
            mem_gamma_direct, mem_gamma_mediated = read_gamma("membrane_gamma.dat")
            m = 1  # membrane environment
            count = 0
            for i in range(20):
                for j in range(i, 20):
                    gamma_ijm[m][i][j] = mem_gamma_direct[count][0]*k_relative_mem
                    gamma_ijm[m][j][i] = mem_gamma_direct[count][0]*k_relative_mem
                    count += 1
            count = 0
            for i in range(20):
                for j in range(i, 20):
                    water_gamma_ijm[m][i][j] = mem_gamma_mediated[count][1]*k_relative_mem
                    water_gamma_ijm[m][j][i] = mem_gamma_mediated[count][1]*k_relative_mem
                    count += 1
            count = 0
            for i in range(20):
                for j in range(i, 20):
                    protein_gamma_ijm[m][i][j] = mem_gamma_mediated[count][0]*k_relative_mem
                    protein_gamma_ijm[m][j][i] = mem_gamma_mediated[count][0]*k_relative_mem
                    count += 1
            for i in range(self.nres):
                for j in range(self.nres):
                    resId1 = i
                    chain1 = inWhichChain(resId1, self.chain_ends)
                    resId2 = j
                    chain2 = inWhichChain(resId2, self.chain_ends)
                    if abs(resId1-resId2)-min_sequence_separation_mem >= 0 or chain1 != chain2:
                        res_table[m][i][j] = 1
                    else:
                        res_table[m][i][j] = 0

        contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, 20, 20, gamma_ijm.T.flatten()))
        contact.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, 20, 20, water_gamma_ijm.T.flatten()))
        contact.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, 20, 20, protein_gamma_ijm.T.flatten()))
        contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(20, 3, burial_gamma.T.flatten()))
        contact.addTabulatedFunction("res_table", Discrete3DFunction(nwell, self.nres, self.nres, res_table.T.flatten()))

        contact.addPerParticleParameter("resName")
        contact.addPerParticleParameter("resId")
        contact.addPerParticleParameter("isCb")
        contact.addGlobalParameter("k_contact", k_contact)
        contact.addGlobalParameter("eta", eta)
        contact.addGlobalParameter("eta_sigma", eta_sigma)
        contact.addGlobalParameter("rho_0", rho_0)
        contact.addGlobalParameter("min_sequence_separation", min_sequence_separation)
        contact.addGlobalParameter("rmin", r_min)
        contact.addGlobalParameter("rmax", r_max)
        contact.addGlobalParameter("rminII", r_minII)
        contact.addGlobalParameter("rmaxII", r_maxII)
        contact.addGlobalParameter("burial_kappa", burial_kappa)

        contact.addComputedValue("rho", "step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)

        # if z_dependent:
        #     contact.addComputedValue("isInMembrane", f"step({z_m}-abs(z))", CustomGBForce.SingleParticle)
        # else:
        #     contact.addComputedValue("isInMembrane", "0", CustomGBForce.SingleParticle)


        # contact.addComputedValue("isInMembrane", "1", CustomGBForce.SingleParticle)
        # replace cb with ca for GLY
        cb_fixed = [x if x > 0 else y for x,y in zip(self.cb,self.ca)]
        none_cb_fixed = [i for i in range(self.natoms) if i not in cb_fixed]
        # print(self.natoms, len(self.resi), self.resi, seq)
        for i in range(self.natoms):
            contact.addParticle([gamma_se_map_1_letter[seq[self.resi[i]]], self.resi[i], int(i in cb_fixed)])


        if z_dependent:
            # print(f"0.5*tanh({eta_switching}*(z+{z_m}))+0.5*tanh({eta_switching}*({z_m}-z))")
            contact.addComputedValue("alphaMembrane", f"0.5*tanh({eta_switching}*(z+{z_m}))+0.5*tanh({eta_switching}*({z_m}-z))", CustomGBForce.SingleParticle)
            # contact.addComputedValue("alphaMembrane", f"z", CustomGBForce.SingleParticle)
            # contact.addComputedValue("isInMembrane", f"z", CustomGBForce.SingleParticle)
            # contact.addComputedValue("isInMembrane", f"step({z_m}-abs(z))", CustomGBForce.SingleParticle)
            # mediated term
            contact.addEnergyTerm("(1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part;\
                                    water_part=-res_table(0, resId1, resId2)*k_contact*thetaII*\
                                    (sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                    sigma_protein*protein_gamma_ijm(0, resName1, resName2));\
                                    membrane_part=-res_table(1, resId1, resId2)*k_contact*thetaII*\
                                    (sigma_water*water_gamma_ijm(1, resName1, resName2)+\
                                    sigma_protein*protein_gamma_ijm(1, resName1, resName2));\
                                    sigma_protein=1-sigma_water;\
                                    thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                    sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                    CustomGBForce.ParticlePair)
            # direct term
            contact.addEnergyTerm("(1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part;\
                                    water_part=-res_table(0, resId1, resId2)*k_contact*\
                                    gamma_ijm(0, resName1, resName2)*theta;\
                                    membrane_part=-res_table(1, resId1, resId2)*k_contact*\
                                    gamma_ijm(1, resName1, resName2)*theta;\
                                    theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
                                    CustomGBForce.ParticlePair)
        else:
            # mediated term
            contact.addEnergyTerm(f"-res_table({inMembrane}, resId1, resId2)*k_contact*thetaII*\
                                    (sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                                    sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2));\
                                    sigma_protein=1-sigma_water;\
                                    thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                    sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                    CustomGBForce.ParticlePair)
            # direct term
            contact.addEnergyTerm(f"-res_table({inMembrane}, resId1, resId2)*k_contact*\
                                    gamma_ijm({inMembrane}, resName1, resName2)*theta;\
                                    theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
                                    CustomGBForce.ParticlePair)

        # burial term
        for i in range(3):
            contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            contact.addEnergyTerm(f"-0.5*isCb*k_contact*burial_gamma_ij(resName, {i})*\
                                        (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                        tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

        print("Number of atom: ", self.natoms, "Number of residue: ", len(cb_fixed))
        # print(len(none_cb_fixed), len(cb_fixed))
        for e1 in none_cb_fixed:
            for e2 in none_cb_fixed:
                if e1 > e2:
                    continue
                contact.addExclusion(e1, e2)
        for e1 in none_cb_fixed:
            for e2 in cb_fixed:
                contact.addExclusion(e1, e2)

        # contact.setCutoffDistance(1.1)
        contact.setNonbondedMethod(CustomGBForce.CutoffNonPeriodic)
        print("Contact cutoff ", contact.getCutoffDistance())
        print("NonbondedMethod: ", contact.getNonbondedMethod())
        contact.setForceGroup(18)
        return contact

    def contact_test_term(self, k_contact=4.184, z_dependent=False, z_m=1.5):
        contact = CustomGBForce()
        gamma_ijm = np.zeros((2, 20, 20))
        contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(2, 20, 20, gamma_ijm.T.flatten()))
        contact.addComputedValue("rho", f"1", CustomGBForce.ParticlePair)
        contact.addComputedValue("alpha", f"z", CustomGBForce.SingleParticle)
        for i in range(self.natoms):
            contact.addParticle()
        return contact


    def fragment_memory_term(self, k_fm=0.04184, frag_location_pre="./",
                        min_seq_sep=3, max_seq_sep=9, fm_well_width=0.1):
        # 0.8368 = 0.01 * 4.184 # in kJ/mol, converted from default value in LAMMPS AWSEM
        k_fm *= self.k_awsem
        frag_table_rmin = 0
        frag_table_rmax = 5 # in nm
        frag_table_dr = 0.01
        r_array = np.arange(frag_table_rmin, frag_table_rmax, frag_table_dr)
        number_of_atoms = self.natoms
        r_table_size = int((frag_table_rmax - frag_table_rmin)/frag_table_dr)  # 500 here.
        raw_frag_table = np.zeros((number_of_atoms, 6*max_seq_sep, r_table_size))
        data_dic = {}
        for i in range(self.natoms):
            if i in self.ca:
                res_id = self.resi[i]    # self.resi start with 0, but pdb residue id start with 1
                data_dic[("CA", 1+int(res_id))] = i
            if i in self.cb:
                res_id = self.resi[i]
                data_dic[("CB", 1+int(res_id))] = i
        # print(self.res_type)
        # print(self.resi)
        # print(data_dic)
        frag_file_list_file = frag_location_pre + "frags.mem"
        frag_table_file = frag_location_pre + "frag_table.npy"

        if os.path.isfile(frag_table_file):
            print(f"Reading Fragment table. from {frag_table_file}.")
            frag_table, interaction_list, interaction_pair_to_bond_index = np.load(frag_table_file)
            print(f"Fragment table loaded, number of bonds: {len(interaction_list)}")
            frag_file_list = []
        else:
            print(f"Fragment table file is not found. Reading fragments files.")
            frag_file_list = pd.read_table(frag_file_list_file, skiprows=4, sep="\s+", names=["location", "target_start", "fragment_start", "frag_len", "weight"])
            interaction_list = set()
        for frag_index in range(len(frag_file_list)):
            location = frag_file_list["location"].iloc[frag_index]
            frag_name = frag_location_pre + location
            frag_len = frag_file_list["frag_len"].iloc[frag_index]
            weight = frag_file_list["weight"].iloc[frag_index]
            target_start = frag_file_list["target_start"].iloc[frag_index] # residue id
            fragment_start = frag_file_list["fragment_start"].iloc[frag_index] # residue id
            frag = pd.read_table(frag_name, skiprows=2, sep="\s+", header=None, names=["Res_id", "Res", "Type", "i", "x", "y", "z"])
            frag = frag.query(f"Res_id >= {fragment_start} and Res_id < {fragment_start+frag_len} and (Type == 'CA' or Type == 'CB')")
            w_m = weight
            gamma_ij = 1
            f = frag.values
            for i in range(len(frag)):
                for j in range(i, len(frag)):
                    res_id_i = frag["Res_id"].iloc[i]
                    res_id_j = frag["Res_id"].iloc[j]
                    target_res_id_i = frag["Res_id"].iloc[i] - fragment_start + target_start
                    target_res_id_j = frag["Res_id"].iloc[j] - fragment_start + target_start
                    seq_sep = res_id_j - res_id_i
                    if seq_sep >= max_seq_sep:
                        continue
                    if seq_sep < min_seq_sep:
                        continue
                    try:
                        i_type = frag["Type"].iloc[i]
                        j_type = frag["Type"].iloc[j]
                        correspond_target_i = data_dic[(i_type, int(target_res_id_i))]
                        correspond_target_j = data_dic[(j_type, int(target_res_id_j))]
                        correspond_target_i = int(correspond_target_i)
                        correspond_target_j = int(correspond_target_j)
                    except Exception as e:
                        continue

                    fi_x = f[i][4]
                    fi_y = f[i][5]
                    fi_z = f[i][6]

                    fj_x = f[j][4]
                    fj_y = f[j][5]
                    fj_z = f[j][6]
                    # print("----", fi_x, fi_y, fi_z, fj_x, fj_y, fj_z)
                    sigma_ij = fm_well_width*seq_sep**0.15
                    rm = ((fi_x-fj_x)**2 + (fi_y-fj_y)**2 + (fi_z-fj_z)**2)**0.5

                    i_j_sep = int(correspond_target_j - correspond_target_i)

                    raw_frag_table[correspond_target_i][i_j_sep] += w_m*gamma_ij*np.exp((r_array-rm)**2/(-2.0*sigma_ij**2))
                    interaction_list.add((correspond_target_i, correspond_target_j))
        if not os.path.isfile(frag_table_file):
            # Reduce memory usage.
            print("Saving fragment table as npy file to speed up future calculation.")
            number_of_bonds = len(interaction_list)
            frag_table = np.zeros((number_of_bonds, r_table_size))
            interaction_pair_to_bond_index = {}
            for index, (i, j) in enumerate(interaction_list):
                ij_sep = j - i
                assert(ij_sep > 0)
                frag_table[index] = raw_frag_table[i][ij_sep]
                interaction_pair_to_bond_index[(i,j)] = index
            np.save(frag_table_file, (frag_table, interaction_list, interaction_pair_to_bond_index))

        # fm = CustomNonbondedForce(f"-k_fm*((v2-v1)*r+v1*r_2-v2*r_1)/(r_2-r_1); \
        #                             v1=frag_table(index_smaller, sep, r_index_1);\
        #                             v2=frag_table(index_smaller, sep, r_index_2);\
        #                             index_smaller=min(index1,index2);\
        #                             sep=abs(index1-index2);\
        #                             r_1=frag_table_rmin+frag_table_dr*r_index_1;\
        #                             r_2=frag_table_rmin+frag_table_dr*r_index_2;\
        #                             r_index_2=r_index_1+1;\
        #                             r_index_1=floor(r/frag_table_dr);")
        # for i in range(self.natoms):
        #     fm.addParticle([i])

        # # add interaction that are cutoff away
        # # print(sorted(interaction_list))
        # for (i, j) in interaction_list:
        #     fm.addInteractionGroup([i], [j])
        # # add per-particle parameters
        # fm.addPerParticleParameter("index")

        fm = CustomCompoundBondForce(2, "-k_fm*((v2-v1)*r+v1*r_2-v2*r_1)/(r_2-r_1); \
                                    v1=frag_table(index, r_index_1);\
                                    v2=frag_table(index, r_index_2);\
                                    r_1=frag_table_rmin+frag_table_dr*r_index_1;\
                                    r_2=frag_table_rmin+frag_table_dr*r_index_2;\
                                    r_index_2=r_index_1+1;\
                                    r_index_1=floor(r/frag_table_dr);\
                                    r=distance(p1, p2);")
        for (i, j) in interaction_list:
            fm.addBond([i, j], [interaction_pair_to_bond_index[(i,j)]])

        fm.addPerBondParameter("index")

        fm.addTabulatedFunction("frag_table",
                Discrete2DFunction(len(interaction_list), r_table_size, frag_table.T.flatten()))
        fm.addGlobalParameter("k_fm", k_fm)
        fm.addGlobalParameter("frag_table_dr", frag_table_dr)
        fm.addGlobalParameter("frag_table_rmin", frag_table_rmin)

        fm.setForceGroup(19)
        return fm


    def membrane_term(self, k_membrane=4.184, k_m=2, z_m=1.5):
        # k_m in units of nm^-1, z_m in units of nm.
        # add membrane forces
        # 1 Kcal = 4.184 kJ strength by overall scaling
        k_membrane *= self.k_awsem
        membrane = CustomExternalForce(f"{k_membrane}*\
                (0.5*tanh({k_m}*(z+{z_m}))+0.5*tanh({k_m}*({z_m}-z)))*hydrophobicityScale")
        membrane.addPerParticleParameter("hydrophobicityScale")
        zim = np.loadtxt("zim")
        cb_fixed = [x if x > 0 else y for x,y in zip(self.cb,self.ca)]
        for i in cb_fixed:
            # print(self.resi[i] , self.seq[self.resi[i]])
            membrane.addParticle(i, [zim[self.resi[i]]])
        membrane.setForceGroup(20)
        return membrane


    def read_memory(self, pdb_file, chain_name, target_start, fragment_start, length, weight, min_seq_sep, max_seq_sep, am_well_width=0.1):
        memory_interactions = []

        # if not os.path.isfile(pdb_file):
        #     pdbl = PDBList()
        #     pdbl.retrieve_pdb_file(pdb_file.split('.')[0].lower(), pdir='.')
        #     os.rename("pdb%s.ent" % pdb_id, "%s.pdb" % pdb_id)

        parser = PDBParser()
        structure = parser.get_structure('X', pdb_file)
        chain = structure[0][chain_name]
        residues = [x for x in chain if x.get_full_id()[3][1] in range(fragment_start,fragment_start+length-1)]
        for i, residue_i in enumerate(residues):
            for j, residue_j in enumerate(residues):
                if abs(i-j) > max_seq_sep:
                    continue
                target_index_i = target_start + i - 1
                target_index_j = target_start + j - 1
                atom_list_i = []
                target_atom_list_i = []
                atom_list_j = []
                target_atom_list_j = []
                if i-j >= min_seq_sep: # taking the signed value to avoid double counting
                    ca_i = residue_i['CA']
                    atom_list_i.append(ca_i)
                    target_atom_list_i.append(self.ca[target_index_i])
                    ca_j = residue_j['CA']
                    atom_list_j.append(ca_j)
                    target_atom_list_j.append(self.ca[target_index_j])
                    if not residue_i.get_resname() == "GLY" and self.cb[target_index_i] >= 0:
                        cb_i = residue_i['CB']
                        atom_list_i.append(cb_i)
                        target_atom_list_i.append(self.cb[target_index_i])
                    if not residue_j.get_resname() == "GLY" and self.cb[target_index_j] >= 0:
                        cb_j = residue_j['CB']
                        atom_list_j.append(cb_j)
                        target_atom_list_j.append(self.cb[target_index_j])
                for atom_i, atom_j in product(atom_list_i, atom_list_j):
                    particle_1 = target_atom_list_i[atom_list_i.index(atom_i)]
                    particle_2 = target_atom_list_j[atom_list_j.index(atom_j)]
                    r_ijm = abs(atom_i - atom_j)/10.0 # convert to nm
                    sigma_ij = am_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    w_m = weight
                    memory_interaction = [particle_1, particle_2, [w_m, gamma_ij, r_ijm, sigma_ij]]
                    memory_interactions.append(memory_interaction)
        return memory_interactions

    def associative_memory_term(self, memories, k_am=0.8368, min_seq_sep=3, max_seq_sep=9, am_well_width=0.1):
        # 0.8368 = 0.2 * 4.184 # in kJ/mol, converted from default value in LAMMPS AWSEM
        #pdbid #chain #target #fragment #length #weight
        # multiply interaction strength by overall scaling
        k_am *= self.k_awsem
        am_function = '-k_am*w_m*gamma_ij*exp(-(r-r_ijm)^2/(2*sigma_ij^2))'
        am = CustomBondForce(am_function)
        am.addGlobalParameter('k_am', k_am)
        am.addPerBondParameter('w_m')
        am.addPerBondParameter('gamma_ij')
        am.addPerBondParameter('r_ijm')
        am.addPerBondParameter('sigma_ij')
        for memory in memories:
            memory_interactions = self.read_memory(*memory, min_seq_sep, max_seq_sep, am_well_width=am_well_width)
            for memory_interaction in memory_interactions:
                am.addBond(*memory_interaction)
        return am



    def density_dependent_associative_memory_term(self, memories, k_am_dd=1.0, am_dd_min_seq_sep=3, am_dd_max_seq_sep=9, eta_density=50, r_density_min=.45, r_density_max=.65, density_alpha=1.0, density_normalization=2.0, rho0=2.6, am_well_width=0.1, density_min_seq_sep=10, density_only_from_native_contacts=False, density_pdb_file=None, density_chain_name=None, density_native_contact_min_seq_sep=4, density_native_contact_threshold=0.8*nanometers):

        k_am_dd *= self.k_awsem

        am_dd = CustomGBForce()

        # add all particles to force
        for i in range(self.natoms):
            am_dd.addParticle([i])

        # add per-particle parameters
        am_dd.addPerParticleParameter("index")

        # add global parameters
        am_dd.addGlobalParameter("k_am_dd", k_am_dd)
        am_dd.addGlobalParameter("eta_density", eta_density)
        am_dd.addGlobalParameter("r_density_min", r_density_min)
        am_dd.addGlobalParameter("r_density_max", r_density_max)
        am_dd.addGlobalParameter("density_alpha", density_alpha)
        am_dd.addGlobalParameter("density_normalization", density_normalization)
        am_dd.addGlobalParameter("rho0", rho0)

        # if density_only_from_native_contacts, read structure to get native contacts
        if density_only_from_native_contacts:
            structure_interactions = self.read_amhgo_structure(pdb_file=density_pdb_file, chain_name=density_chain_name, amhgo_min_seq_sep=density_native_contact_min_seq_sep, amhgo_contact_threshold=density_native_contact_threshold, amhgo_well_width=0.1) # the well width is not used, so the value doesn't matter

            native_contacts = []
            for interaction in structure_interactions:
                i_index, j_index, [gamma_ij, r_ijN, sigma_ij] = interaction
                native_contacts.append((i_index, j_index))
                native_contacts.append((j_index, i_index))

        # setup tabulated functions and interactions
        density_gamma_ij = [0.0]*self.natoms*self.natoms
        for i in range(self.natoms):
            for j in range(self.natoms):
                if (i in self.cb or (self.res_type[self.resi[i]] == "IGL" and i in self.ca)) and (j in self.cb or (self.res_type[self.resi[j]] == "IGL" and i in self.ca)) and abs(self.resi[i]-self.resi[j])>=density_min_seq_sep:
                    if not density_only_from_native_contacts or (i, j) in native_contacts or (j, i) in native_contacts:
                        density_gamma_ij[i+j*self.natoms] = 1.0
                        density_gamma_ij[j+i*self.natoms] = 1.0
        am_dd.addTabulatedFunction("density_gamma_ij", Discrete2DFunction(self.natoms, self.natoms, density_gamma_ij))

        gamma_ij = [0.0]*self.natoms*self.natoms*len(memories)
        sigma_ij = [0.1]*self.natoms*self.natoms*len(memories)
        r_ijm = [0.0]*self.natoms*self.natoms*len(memories)
        for k, memory in enumerate(memories):
            memory_interactions = self.read_memory(*memory, am_dd_min_seq_sep, am_dd_max_seq_sep, am_well_width=am_well_width)
            for memory_interaction in memory_interactions:
                i, j, (w_m, gamma, r, sigma) = memory_interaction
                gamma_ij[i+j*self.natoms+k*self.natoms*self.natoms] = gamma
                gamma_ij[j+i*self.natoms+k*self.natoms*self.natoms] = gamma
                sigma_ij[i+j*self.natoms+k*self.natoms*self.natoms] = sigma
                sigma_ij[j+i*self.natoms+k*self.natoms*self.natoms] = sigma
                r_ijm[i+j*self.natoms+k*self.natoms*self.natoms] = r
                r_ijm[j+i*self.natoms+k*self.natoms*self.natoms] = r
        am_dd.addTabulatedFunction("gamma_ij", Discrete3DFunction(self.natoms, self.natoms, len(memories), gamma_ij))
        am_dd.addTabulatedFunction("sigma_ij", Discrete3DFunction(self.natoms, self.natoms, len(memories), sigma_ij))
        am_dd.addTabulatedFunction("r_ijm", Discrete3DFunction(self.natoms, self.natoms, len(memories), r_ijm))

        # add computed values
        # compute the density
        am_dd.addComputedValue("rho", "0.25*density_gamma_ij(index1, index2)*(1+tanh(eta_density*(r-r_density_min)))*(1+tanh(eta_density*(r_density_max-r)))", CustomGBForce.ParticlePair)

        # function that determines how the AM term depends on density
        #f_string = "0.25*(1-tanh(eta_density*(rho0-rho1)))*(1-tanh(eta_density*(rho0-rho2)))" # both residues must be buried for the interaction to be active
        f_string = "1-(0.25*(1-tanh(eta_density*(rho1-rho0)))*(1-tanh(eta_density*(rho2-rho0))))" # one residue being buried is enough for the interaction to be active

        # add energy term for each memory
        for k, memory in enumerate(memories):
            memory_interactions = self.read_memory(*memory, am_dd_min_seq_sep, am_dd_max_seq_sep, am_well_width=am_well_width)
            for memory_interaction in memory_interactions:
                i, j, (w_m, gamma, r, sigma) = memory_interaction
            am_dd.addEnergyTerm("-k_am_dd*(density_alpha*f*density_normalization*beta_ij+(1-density_alpha)*beta_ij);\
            beta_ij=%f*gamma_ij(index1,index2,%d)*exp(-(r-r_ijm(index1,index2,%d))^2/(2*sigma_ij(index1,index2,%d)^2));\
            f=%s" % (w_m, k, k, k, f_string), CustomGBForce.ParticlePair)

        return am_dd

    def read_beta_parameters(self):
        ### directly copied from Nick Schafer's
        #os.chdir(parameter_directory)
        in_anti_HB = open("anti_HB", 'r').readlines()
        in_anti_NHB = open("anti_NHB", 'r').readlines()
        in_para_HB = open("para_HB", 'r').readlines()
        in_para_one = open("para_one", 'r').readlines()
        in_anti_one = open("anti_one", 'r').readlines()

        p_par = np.zeros((20))
        p_anti = np.zeros((20))
        p_antihb = np.zeros((20,20,2))
        p_antinhb = np.zeros((20,20,2))
        p_parhb = np.zeros((20,20,2))

        for i in range(20):
            p_par[i] = float(in_para_one[i].strip())
            p_anti[i] = float(in_anti_one[i].strip())
            for j in range(20):
                p_antihb[i][j][0] = float(in_anti_HB[i].strip().split()[j])
                p_antinhb[i][j][0] = float(in_anti_NHB[i].strip().split()[j])
                p_parhb[i][j][0] = float(in_para_HB[i].strip().split()[j])

        for i in range(20):
            for j in range(20):
                p_antihb[i][j][1] = float(in_anti_HB[i+21].strip().split()[j])
                p_antinhb[i][j][1] = float(in_anti_NHB[i+21].strip().split()[j])
                p_parhb[i][j][1] = float(in_para_HB[i+21].strip().split()[j])
        return p_par, p_anti, p_antihb, p_antinhb, p_parhb

    def lambda_coefficient(self, i, j, lambda_index):
        p_par, p_anti, p_antihb, p_antinhb, p_parhb = self.read_beta_parameters()
        parameter_i = []
        #print(i,j,lambda_index)
        for ii in range(self.nres):
            #print(self.seq[i])
            parameter_i.append(se_map_1_letter[self.seq[ii]])
        #print(p_antihb[parameter_i[i], parameter_i[j]][0],p_antinhb[parameter_i[i+1],parameter_i[j-1]][0],p_anti[parameter_i[i]], p_anti[parameter_i[j]])
        lambda_2_extra_terms = -0.5*self.alpha_coefficient(parameter_i[i],parameter_i[j],1)*p_antihb[parameter_i[i], parameter_i[j]][0]-0.25*self.alpha_coefficient(parameter_i[i], parameter_i[j], 2)*(p_antinhb[parameter_i[i+1],parameter_i[j-1]][0] + p_antinhb[parameter_i[i-1],parameter_i[j+1]][0])-self.alpha_coefficient(parameter_i[i], parameter_i[j], 3)*(p_anti[parameter_i[i]]+p_anti[parameter_i[j]])
        lambda_3_extra_terms = -self.alpha_coefficient(parameter_i[i],parameter_i[j], 4)*p_parhb[parameter_i[i+1],parameter_i[j]][0]-self.alpha_coefficient(parameter_i[i],parameter_i[j],5)*p_par[parameter_i[i+1]]+self.alpha_coefficient(parameter_i[i],parameter_i[j],4)*p_par[parameter_i[j]]
        if abs(j-i) >= 4 and abs(j-i) < 18:
            if lambda_index == 1:
                return 1.37
            elif lambda_index == 2:
                return 3.89+lambda_2_extra_terms
            elif lambda_index == 3:
                return 0.0+lambda_3_extra_terms
        elif abs(j-i) >= 18 and abs(j-i) < 45:
            if lambda_index == 1:
                return 1.36
            elif lambda_index == 2:
                return 3.50+lambda_2_extra_terms
            elif lambda_index == 3:
                return 3.47+lambda_3_extra_terms
        elif abs(j-i) >= 45:
            if lambda_index == 1:
                return 1.17
            elif lambda_index == 2:
                return 3.52+lambda_2_extra_terms
            elif lambda_index == 3:
                return 3.62+lambda_3_extra_terms
        elif abs(j-i) < 4:
            return 0.0

    def alpha_coefficient(self, i,j, alpha_index):
        if abs(j-i) >= 4 and abs(j-i) < 18:
            if alpha_index == 1:
                return 1.3
            if alpha_index == 2:
                return 1.32
            if alpha_index == 3:
                return 1.22
            if alpha_index == 4:
                return 0.0
            if alpha_index == 5:
                return 0.0
        elif abs(j-i) >= 18 and abs(j-i) < 45:
            if alpha_index == 1:
                return 1.3
            if alpha_index == 2:
                return 1.32
            if alpha_index == 3:
                return 1.22
            if alpha_index == 4:
                return 0.33
            if alpha_index == 5:
                return 1.01
        elif abs(j-i) >= 45:
            if alpha_index == 1:
                return 1.3
            if alpha_index == 2:
                return 1.32
            if alpha_index == 3:
                return 1.22
            if alpha_index == 4:
                return 0.33
            if alpha_index == 5:
                return 1.01
        elif abs(j-i) <4:
            return 0.0

    def get_lambda_by_index(self, i, j, lambda_i):


        lambda_table = [[1.37, 1.36, 1.17],
                        [3.89, 3.50, 3.52],
                        [0.0,  3.47, 3.62]]
        if abs(j-i) >= 4 and abs(j-i) < 18:
            return lambda_table[lambda_i][0]
        elif abs(j-i) >= 18 and abs(j-i) < 45:
            return lambda_table[lambda_i][1]
        elif abs(j-i) >= 45:
            return lambda_table[lambda_i][2]
        else:
            return 0

    def get_alpha_by_index(self, i, j, alpha_i):
        alpha_table = [[1.30, 1.30, 1.30],
                        [1.32, 1.32, 1.32],
                        [1.22, 1.22, 1.22],
                        [0,    0.33, 0.33],
                        [0.0,  1.01, 1.01]]
        if abs(j-i) >= 4 and abs(j-i) < 18:
            return alpha_table[alpha_i][0]
        elif abs(j-i) >= 18 and abs(j-i) < 45:
            return alpha_table[alpha_i][1]
        elif abs(j-i) >= 45:
            return alpha_table[alpha_i][2]
        else:
            return 0

    def get_Lambda_2(self, i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb):
        Lambda = self.get_lambda_by_index(i, j, 1)
        a = []
        for ii in range(self.nres):
            a.append(se_map_1_letter[self.seq[ii]])

        Lambda += -0.5*self.get_alpha_by_index(i, j, 0)*p_antihb[a[i], a[j]][0]
        Lambda += -0.25*self.get_alpha_by_index(i, j, 1)*(p_antinhb[a[i+1], a[j-1]][0] + p_antinhb[a[i-1], a[j+1]][0])
        Lambda += -self.get_alpha_by_index(i, j, 2)*(p_anti[a[i]] + p_anti[a[j]])
        return Lambda

    def get_Lambda_3(self, i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb):
        Lambda = self.get_lambda_by_index(i, j, 2)
        a = []
        for ii in range(self.nres):
            a.append(se_map_1_letter[self.seq[ii]])

        Lambda += -self.get_alpha_by_index(i, j, 3)*p_parhb[a[i+1], a[j]][0]
        Lambda += -self.get_alpha_by_index(i, j, 4)*p_par[a[i+1]]
        Lambda += -self.get_alpha_by_index(i, j, 3)*p_par[a[j]]
        return Lambda

    def apply_beta_term_1(self, k_beta=4.184):

        print("beta_1 term ON");
        nres, n, h, ca, o, res_type = self.nres, self.n, self.h, self.ca, self.o, self.res_type
        #print(lambda_1)
        r_ON = .298
        sigma_NO = .068
        r_OH = .206
        sigma_HO = .076
        eta_beta_1 = 10.0
        eta_beta_2 = 5.0
        # r_HB_c = 0.4
        r_HB_c = 1.2

        theta_ij =   f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
        # theta_ji =   f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hi-{r_OH})^2/(2*{sigma_HO}^2))"
        # theta_jip2 = "exp(-(r_Oj_Nip2-r_ON)^2/(2*sigma_NO^2)-(r_Oj_Hip2-r_OH)^2/(2*sigma_HO^2))"
        nu_i = f"0.5*(1+tanh({eta_beta_1}*(r_CAim2_CAip2-{r_HB_c})))"
        nu_j = f"0.5*(1+tanh({eta_beta_2}*(r_CAjm2_CAjp2-{r_HB_c})))"

        # Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4     5     6     7
        beta_string_1 = f"-k_beta*lambda_1*theta_ij*nu_i*nu_j;theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                        nu_i={nu_i};nu_j={nu_j};r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)"
        # beta_string_1 = f"-k_beta*lambda_1"
        # beta_string_1 = f"-k_beta"

        beta_1 = CustomCompoundBondForce(7, beta_string_1)
        #beta_2 = CustomCompoundBondForce(10, beta_string_2)
        #beta_3 = CustomCompoundBondForce(10, beta_string_3)
        # add parameters to force
        beta_1.addGlobalParameter("k_beta", k_beta)
        beta_1.addPerBondParameter("lambda_1")
        #beta_2.addTabulatedFunction("lambda_2", Discrete2DFunction(nres, nres, lambda_2))
        #beta_3.addTabulatedFunction("lambda_3", Discrete2DFunction(nres, nres, lambda_3))

        for i in range(nres):
            for j in range(nres):
                if isChainEdge(i, self.chain_starts, self.chain_ends, n=2) or \
                    isChainEdge(j, self.chain_starts, self.chain_ends, n=2):
                    continue
                if not res_type[j] == "IPR":
                    beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [self.get_lambda_by_index(i, j, 0)])
                #if not res_type[i] == "IPR" and not res_type[j] == "IPR":
                #    beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
                #if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
                #    beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])

        beta_1.setForceGroup(23)
        #beta_2.setForceGroup(24)
        #beta_3.setForceGroup(25)
        return beta_1

    def apply_beta_term_2(self, k_beta=4.184):
        print("beta_2 term ON");
        nres, n, h, ca, o, res_type = self.nres, self.n, self.h, self.ca, self.o, self.res_type
        # add beta potential
        # setup parameters
        r_ON = .298
        sigma_NO = .068
        r_OH = .206
        sigma_HO = .076
        eta_beta_1 = 10.0
        eta_beta_2 = 5.0
        # r_HB_c = 0.4
        r_HB_c = 1.2
        p_par, p_anti, p_antihb, p_antinhb, p_parhb = self.read_beta_parameters()

        theta_ij =   f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
        theta_ji =   f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hi-{r_OH})^2/(2*{sigma_HO}^2))"
        nu_i = f"0.5*(1+tanh({eta_beta_1}*(r_CAim2_CAip2-{r_HB_c})))"
        nu_j = f"0.5*(1+tanh({eta_beta_2}*(r_CAjm2_CAjp2-{r_HB_c})))"

        # Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4     5     6     7
        #beta_string_1 = "-k_beta*lambda_1(index_i,index_j)*theta_ij*nu_i*nu_j;theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
        #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)" % (theta_ij, nu_i, nu_j)

        # Oi Nj Hj Oj Ni Hi CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4  5  6  7     8     9     10
        beta_string_2 = f"-k_beta*lambda_2*theta_ij*theta_ji*nu_i*nu_j;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                        theta_ji={theta_ji};r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
                        nu_i={nu_i};nu_j={nu_j};r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)"

        # Oi Nj Hj Oj Ni+2 Hi+2 CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4  5    6    7     8     9     10
        #beta_string_3 = "-k_beta*lambda_3(index_i,index_j)*theta_ij*theta_jip2*nu_i*nu_j;\
        #                theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
        #                theta_ji=%s;r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
        #                theta_jip2=%s;r_Oj_Nip2=distance(p4,p5);r_Oj_Hip2=distance(p4,p6);\
        #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)" % (theta_ij, theta_ji, theta_jip2, nu_i, nu_j)

        #beta_1 = CustomCompoundBondForce(7, beta_string_1)
        beta_2 = CustomCompoundBondForce(10, beta_string_2)
        #beta_3 = CustomCompoundBondForce(10, beta_string_3)
        # add parameters to force
        beta_2.addGlobalParameter("k_beta", k_beta)
        beta_2.addPerBondParameter("lambda_2")



        for i in range(nres):
            for j in range(nres):
                if isChainEdge(i, self.chain_starts, self.chain_ends, n=2) or \
                    isChainEdge(j, self.chain_starts, self.chain_ends, n=2):
                    continue
                #if not res_type[j] == "IPR":
                #    beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
                if not res_type[i] == "IPR" and not res_type[j] == "IPR":
                    beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [self.get_Lambda_2(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb)])
                #if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
                #    beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])


        #beta_1.setForceGroup(23)
        beta_2.setForceGroup(24)
        #beta_3.setForceGroup(25)
        return beta_2

    def apply_beta_term_3(self, k_beta=4.184):
        print("beta_3 term ON");
        nres, n, h, ca, o, res_type = self.nres, self.n, self.h, self.ca, self.o, self.res_type
        # add beta potential
        # setup parameters
        r_ON = .298
        sigma_NO = .068
        r_OH = .206
        sigma_HO = .076
        eta_beta_1 = 10.0
        eta_beta_2 = 5.0
        # r_HB_c = 0.4
        r_HB_c = 1.2
        p_par, p_anti, p_antihb, p_antinhb, p_parhb = self.read_beta_parameters()

        theta_ij =   f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
        theta_jip2 = f"exp(-(r_Oj_Nip2-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hip2-{r_OH})^2/(2*{sigma_HO}^2))"
        nu_i = f"0.5*(1+tanh({eta_beta_1}*(r_CAim2_CAip2-{r_HB_c})))"
        nu_j = f"0.5*(1+tanh({eta_beta_2}*(r_CAjm2_CAjp2-{r_HB_c})))"

        # Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4     5     6     7
        #beta_string_1 = "-k_beta*lambda_1(index_i,index_j)*theta_ij*nu_i*nu_j;theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
        #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)" % (theta_ij, nu_i, nu_j)

        # Oi Nj Hj Oj Ni Hi CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4  5  6  7     8     9     10
        #beta_string_2 = "-k_beta*lambda_2(index_i,index_j)*theta_ij*theta_ji*nu_i*nu_j;\
        #                theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
        #                theta_ji=%s;r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
        #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)" % (theta_ij, theta_ji, nu_i, nu_j)

        # Oi Nj Hj Oj Ni+2 Hi+2 CAi-2 CAi+2 CAj-2 CAj+2
        # 1  2  3  4  5    6    7     8     9     10
        beta_string_3 = f"-k_beta*lambda_3*theta_ij*theta_jip2*nu_i*nu_j;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                        theta_jip2={theta_jip2};r_Oj_Nip2=distance(p4,p5);r_Oj_Hip2=distance(p4,p6);\
                        nu_i={nu_i};nu_j={nu_j};r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)"

        beta_3 = CustomCompoundBondForce(10, beta_string_3)
        # add parameters to force
        beta_3.addGlobalParameter("k_beta", k_beta)
        beta_3.addPerBondParameter("lambda_3")

        for i in range(nres):
            for j in range(nres):
                if isChainEdge(i, self.chain_starts, self.chain_ends, n=2) or \
                    isChainEdge(j, self.chain_starts, self.chain_ends, n=2):
                    continue
                #if not res_type[j] == "IPR":
                #    beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
                #if not res_type[i] == "IPR" and not res_type[j] == "IPR":
                #    beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
                if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
                    beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [self.get_Lambda_3(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb)])


        #beta_1.setForceGroup(23)
        #beta_2.setForceGroup(24)
        beta_3.setForceGroup(25)
        return beta_3

    def pap_term(self, k_pap=4.184):
        print("pap term ON");
        nres, ca = self.nres, self.ca
        # r0 = 2.0 # nm
        r0 = 0.8 # nm
        eta_pap = 70 # nm^-1
        gamma_aph = 1.0
        gamma_ap = 0.4
        gamma_p = 0.4
        pap_function = f"-k_pap*gamma*0.5*(1+tanh({eta_pap}*({r0}-distance(p1,p2))))*0.5*(1+tanh({eta_pap}*({r0}-distance(p3,p4))))"
        pap = CustomCompoundBondForce(4, pap_function)
        pap.addGlobalParameter("k_pap", k_pap)
        pap.addPerBondParameter("gamma")
        #count = 0;
        for i in range(nres):
            for j in range(nres):
                # anti-parallel hairpin for i from 1 to N-13 and j from i+13 to min(i+16,N)
                # CAi CAj CAi+4 CAj-4
                # 1   2   3     4
                if i <= nres-13 and j >= i+13 and j <= min(i+16,nres):
                    pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_aph])
                    #count = count + 1
                    #print([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_aph])
                # anti-parallel for i from 1 to N-17 and j from i+17 to N
                # CAi CAj CAi+4 CAj-4
                # 1   2   3     4
                if i <= nres-17 and j >= i+17 and j <= nres:
                    pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_ap])
                    #count = count + 1;
                    #print([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_ap])
                # parallel for i from 1 to N-13 and j from i+9 to N-4
                # CAi CAj CAi+4 CAj+4
                # 1   2   3     4
                if i <= nres-13 and j >= i+9 and j < nres-4:
                    #print([i, j, i+4, j+4])
                    #print([i, j, i+4, j+4, ca[i], ca[j], ca[i+4], ca[j+4]], [gamma_p])
                    pap.addBond([ca[i], ca[j], ca[i+4], ca[j+4]], [gamma_p])
                    #count = count + 1;

        #print(count)
        pap.setForceGroup(26)
        return pap

    def read_amhgo_structure(self, pdb_file, chain_name, amhgo_min_seq_sep=4, amhgo_contact_threshold=0.8*nanometers, amhgo_well_width=0.1):
        structure_interactions = []
        parser = PDBParser()
        structure = parser.get_structure('X', pdb_file)
        chain = structure[0][chain_name]
        residues = [x for x in chain]
        for i, residue_i in enumerate(residues):
            for j, residue_j in enumerate(residues):
                ca_list = []
                cb_list = []
                atom_list_i = []
                atom_list_j = []
                if i-j >= amhgo_min_seq_sep:  # taking the signed value to avoid double counting
                    ca_i = residue_i['CA']
                    ca_list.append(ca_i)
                    atom_list_i.append(ca_i)
                    ca_j = residue_j['CA']
                    ca_list.append(ca_j)
                    atom_list_j.append(ca_j)
                    if not residue_i.get_resname() == "GLY":
                        cb_i = residue_i['CB']
                        cb_list.append(cb_i)
                        atom_list_i.append(cb_i)
                    if not residue_j.get_resname() == "GLY":
                        cb_j = residue_j['CB']
                        cb_list.append(cb_j)
                        atom_list_j.append(cb_j)
                    for atom_i, atom_j in product(atom_list_i, atom_list_j):
                        r_ijN = abs(atom_i - atom_j)/10.0*nanometers # convert to nm
                        if r_ijN <= amhgo_contact_threshold:
                            sigma_ij = amhgo_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                            gamma_ij = 1.0
                            if atom_i in ca_list:
                                i_index = self.ca[i]
                            if atom_i in cb_list:
                                i_index = self.cb[i]
                            if atom_j in ca_list:
                                j_index = self.ca[j]
                            if atom_j in cb_list:
                                j_index = self.cb[j]
                            structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                            print(i_index, j_index, gamma_ij, r_ijN, sigma_ij)
                            structure_interactions.append(structure_interaction)
        return structure_interactions

    def additive_amhgo_term(self, pdb_file, chain_name, k_amhgo=4.184, amhgo_min_seq_sep=10, amhgo_contact_threshold=0.8*nanometers, amhgo_well_width=0.1):
        import itertools
        # multiply interaction strength by overall scaling
        print("AMH-GO structure based term is ON")
        k_amhgo *= self.k_awsem
        # create contact force
        amhgo = CustomBondForce("-k_amhgo*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
        # # add global parameters
        amhgo.addGlobalParameter("k_amhgo", k_amhgo)
        amhgo.addPerBondParameter("gamma_ij")
        amhgo.addPerBondParameter("r_ijN")
        amhgo.addPerBondParameter("sigma_ij")
        # create bonds
        structure_interactions = self.read_amhgo_structure(pdb_file, chain_name, amhgo_min_seq_sep, amhgo_contact_threshold, amhgo_well_width=amhgo_well_width)
        print(structure_interactions)
        for structure_interaction in structure_interactions:
            print(structure_interaction)
            amhgo.addBond(*structure_interaction)
        #amhgo.setForceGroup(22)
        return amhgo

    def er_term(self, k_er=4.184, er_min_seq_sep=2, er_cutoff=99.0, er_well_width=0.1):
        ### this is a structure prediction related term; Adapted from Sirovitz Schafer Wolynes 2017 Protein Science;
        ### See original papers for reference: Make AWSEM AWSEM-ER with Evolutionary restrictions
        ### ER restrictions can be obtained from multiple sources (RaptorX, deepcontact, and Gremlin)
        ### term modified from amh-go term, and the current strength seems to be high, and needs to be lowered somehow.
        ### amh-go normalization factor will be added soon. Based on Eastwood Wolynes 2000 JCP
        print("ER term is ON")
        import itertools
        k_er *= self.k_awsem
        # create contact force
        er = CustomBondForce("-k_er*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
        # # add global parameters
        er.addGlobalParameter("k_er", k_er)
        er.addPerBondParameter("gamma_ij")
        er.addPerBondParameter("r_ijN")
        er.addPerBondParameter("sigma_ij")
        structure_interactions_er = []
        ### read in dat files from contact predictions;
        in_rnativeCACA = np.loadtxt('go_rnativeCACA.dat');
        in_rnativeCACB = np.loadtxt('go_rnativeCACB.dat');
        in_rnativeCBCB = np.loadtxt('go_rnativeCBCB.dat');
        for i in range(self.nres):
            for j in range(self.nres):
                if abs(i-j) >= er_min_seq_sep and in_rnativeCACA[i][j]<er_cutoff:
                    sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    r_ijN = in_rnativeCACA[i][j]/10.0*nanometers;
                    structure_interactions_er.append([self.ca[i], self.ca[j], [gamma_ij, r_ijN, sigma_ij]])
                if abs(i-j) >= er_min_seq_sep and in_rnativeCACB[i][j]<er_cutoff and self.cb[j]!= -1:
                    sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    r_ijN = in_rnativeCACB[i][j]/10.0*nanometers;
                    structure_interactions_er.append([self.ca[i], self.cb[j], [gamma_ij, r_ijN, sigma_ij]])
                if abs(i-j) >= er_min_seq_sep and in_rnativeCBCB[i][j]<er_cutoff and self.cb[j]!= -1 and self.cb[i]!= -1:#self.res_type[self.resi[i]] != "IGL" and self.res_type[self.resi[j]] != "IGL":
                    sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    r_ijN = in_rnativeCBCB[i][j]/10.0*nanometers;
                    structure_interactions_er.append([self.cb[i], self.cb[j], [gamma_ij, r_ijN, sigma_ij]])
                    #print([i, j, self.res_type[self.resi[i]], self.res_type[self.resi[j]],self.cb[i], self.cb[j], [gamma_ij, r_ijN, sigma_ij]])
        # create bonds
        for structure_interaction_er in structure_interactions_er:
            er.addBond(*structure_interaction_er)
        er.setForceGroup(21)
        return er

    def debye_huckel_term(self, k_dh = 4.15*4.184):
         print("Debye Huckel term is ON")
         k_dh *= self.k_awsem*0.1;
         k_screening = 1.0;
         screening_length = 1.0; #(in the unit of nanometers)

         dh = CustomBondForce("k_dh*charge_i*charge_j/r*exp(-k_screening*r/screening_length)")
         dh.addGlobalParameter("k_dh", k_dh)
         dh.addGlobalParameter("k_screening", k_screening)
         dh.addGlobalParameter("screening_length", screening_length)
         dh.addPerBondParameter("charge_i");
         dh.addPerBondParameter("charge_j")
         structure_interactions_dh = []
         for i in range(self.nres):
             for j in range(i+1,self.nres):
                 charge_i = 0.0;
                 charge_j = 0.0;
                 if self.seq[i] == "R" or self.seq[i]=="K":
                     charge_i = 1.0;
                 if self.seq[i] == "D" or self.seq[i]=="E":
                     charge_i = -1.0;
                 if self.seq[j] == "R" or self.seq[j]=="K":
                     charge_j = 1.0;
                 if self.seq[j] == "D" or self.seq[j]=="E":
                     charge_j = -1.0;
                 if charge_i*charge_j!=0.0:
                     structure_interactions_dh.append([self.cb[i], self.cb[j], [charge_i, charge_j]]);
                     #print([self.seq[i], self.seq[j],self.cb[i], self.cb[j], [charge_i, charge_j]])
         for structure_interaction_dh in structure_interactions_dh:
             dh.addBond(*structure_interaction_dh)
         dh.setForceGroup(27)
         return dh

    def qbias_term(self, q0, reference_pdb_file, reference_chain_name, k_qbias=10000, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.8*nanometers):
        qbias = CustomCVForce("0.5*k_qbias*(q-q0)^2")
        q = self.q_value(reference_pdb_file, reference_chain_name, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
        qbias.addCollectiveVariable("q", q)
        qbias.addGlobalParameter("k_qbias", k_qbias)
        qbias.addGlobalParameter("q0", q0)
        return qbias

def read_trajectory_pdb_positions(pdb_trajectory_filename):
    import uuid, os
    pdb_trajectory_contents = open(pdb_trajectory_filename).read().split("MODEL")[1:]
    pdb_trajectory_contents = ['\n'.join(x.split('\n')[1:]) for x in pdb_trajectory_contents]
    pdb_trajectory = []
    for i, pdb_contents in enumerate(pdb_trajectory_contents):
        temporary_file_name = str(uuid.uuid4())
        temporary_pdb_file = open(temporary_file_name, 'w')
        temporary_pdb_file.write(pdb_contents)
        temporary_pdb_file.close()
        pdb = PDBFile(temporary_file_name)
        pdb_trajectory.append(pdb)
        os.remove(temporary_file_name)
    return pdb_trajectory

def compute_order_parameters(openmm_awsem_pdb_file, pdb_trajectory_filename, order_parameters,
            platform_name='CPU', k_awsem=1.0, compute_mdtraj=False, rmsd_reference_structure=None,
            compute_total_energy=True, energy_columns=None, xml_filename="awsem.xml"):
    pdb_trajectory = read_trajectory_pdb_positions(pdb_trajectory_filename)
    order_parameter_values = []
    for i, order_parameter in enumerate(order_parameters):
        order_parameter_values.append([])
        oa = OpenMMAWSEMSystem(openmm_awsem_pdb_file, k_awsem=k_awsem, xml_filename=xml_filename)
        platform = Platform.getPlatformByName(platform_name) # OpenCL, CUDA, CPU, or Reference
        integrator = VerletIntegrator(2*femtoseconds)
        oa.addForce(order_parameter)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        for pdb in pdb_trajectory:
            simulation.context.setPositions(pdb.positions)
            state = simulation.context.getState(getEnergy=True)
            order_parameter_values[i].append(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
    if compute_mdtraj:
        md_traj_order_parameters = compute_mdtraj_order_parmeters(pdb_trajectory_filename, rmsd_reference_structure=rmsd_reference_structure)
        for key, value in md_traj_order_parameters.items():
            if len(value.flatten()) == len(pdb_trajectory):
                order_parameter_values.append(value)
    if compute_total_energy:
        order_parameter_values.append([])
        for i in range(len(pdb_trajectory)):
            order_parameter_values[-1].append(np.sum([order_parameter_values[x-1][i] for x in energy_columns]))
    if not compute_mdtraj:
        return np.array(order_parameter_values)
    else:
        return np.array(order_parameter_values), md_traj_order_parameters

def compute_mdtraj_order_parmeters(trajectory_file, rmsd_reference_structure=None):
    # documentation: http://mdtraj.org/1.8.0/analysis.html#
    trajectory = md.load(trajectory_file)

    return_values = []
    return_value_names = []

    if not rmsd_reference_structure == None:
        reference = md.load(rmsd_reference_structure)
        rmsd = md.rmsd(trajectory, reference)
        return_values.append(rmsd)
        return_value_names.append("RMSD")

    hydrogen_bonds = np.array([np.sum(x) for x in md.kabsch_sander(trajectory)])
    return_values.append(hydrogen_bonds)
    return_value_names.append("HBondEnergy")

    ss = md.compute_dssp(trajectory)
    shape = ss.shape
    transdict = dict(zip(list(set(list(ss.flatten()))),range(len(list(set(list(ss.flatten())))))))
    ss = np.array([transdict[x] for x in ss.flatten()]).reshape(shape).T
    return_values.append(ss)
    return_value_names.append("SecondaryStructure")

    rg = md.compute_rg(trajectory)
    return_values.append(rg)
    return_value_names.append("Rg")

    distances, residue_pairs = md.compute_contacts(trajectory, scheme='ca')
    contacts = md.geometry.squareform(distances, residue_pairs)
    return_values.append(contacts)
    return_value_names.append("Contacts")

    return dict(zip(return_value_names, return_values))

def plot_free_energy(order_parameter_file, labels, bins=20, discard=.2, twodlimits=[[0,1],[0,1]]):
    import matplotlib as mpl
    import scipy.ndimage
    mpl.rcParams['font.size'] = 24

    if len(labels) > 2:
        print("Too many labels to plot.")
        return
    order_parameters = np.loadtxt(order_parameter_file).T
    data_labels = open(order_parameter_file).readlines()[0].split()[1:]
    data = []
    for label in labels:
        raw_data = order_parameters[data_labels.index(label)]
        data_after_discard = raw_data[int(len(raw_data)*discard):]
        comparison_data_set_size = (1.0-discard)/2
        data_set_1 = raw_data[int(len(raw_data)*discard):int(len(raw_data)*(discard+comparison_data_set_size))]
        data_set_2 = raw_data[int(len(raw_data)*(discard+comparison_data_set_size)):]
        data.append(data_after_discard)
    if len(labels) == 1:
        hist, bins = np.histogram(data[0], density=True, bins=bins)
        hist /= np.sum(hist)
        f = -np.log(hist)
        bin_centers = (bins[:-1]+bins[1:])/2
        hist, bins = np.histogram(data_set_1, density=True, bins=bins)
        hist /= np.sum(hist)
        f_1 = -np.log(hist)
        bin_centers_1 = (bins[:-1]+bins[1:])/2
        hist, bins = np.histogram(data_set_2, density=True, bins=bins)
        hist /= np.sum(hist)
        f_2 = -np.log(hist)
        bin_centers_2 = (bins[:-1]+bins[1:])/2
        plt.figure(figsize=(12,8))
        plt.plot(bin_centers, f, label="All data")
        plt.plot(bin_centers_1, f_1, label="Split data 1")
        plt.plot(bin_centers_2, f_2, label="Split data 2")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title(labels[0])
        plt.xlabel(labels[0])
        plt.ylabel("Free Energy (kT)")
    if len(labels) == 2:
        H, xedges, yedges = np.histogram2d(data[0], data[1], bins=bins)
        xcenters = (xedges[:-1] + xedges[1:])/2.0
        ycenters = (yedges[:-1] + yedges[1:])/2.0
        H /= np.sum(H)
        H = -np.log(H)
        H -= np.min(H)

        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111)
        ax.set_title("%s vs. %s" % (labels[0], labels[1]))
        ax.set_xlim(twodlimits[0])
        ax.set_ylim(twodlimits[1])
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        X, Y = np.meshgrid(xcenters, ycenters)
        plt.contourf(X, Y, H, cmap=plt.cm.get_cmap('jet'))

        plt.colorbar()
        plt.tight_layout()
        plt.show()

def compute_perturbed_energies(openmm_awsem_pdb_file, pdb_trajectory_filename, perturbations, order_parameter_values, platform_name='CPU', k_awsem=1.0, total_energy_column=15):
    pdb_trajectory = read_trajectory_pdb_positions(pdb_trajectory_filename)
    all_perturbed_energies = []
    for i, perturbation in enumerate(perturbations):
        # compute new energies
        perturbed_energies = []
        for j, energy_term in enumerate(perturbation):
            perturbed_energies.append([])
            oa = OpenMMAWSEMSystem(openmm_awsem_pdb_file, k_awsem=k_awsem)
            platform = Platform.getPlatformByName(platform_name) # OpenCL, CUDA, CPU, or Reference
            integrator = VerletIntegrator(2*femtoseconds)
            oa.addForce(energy_term[1])
            simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
            for pdb in pdb_trajectory:
                simulation.context.setPositions(pdb.positions)
                state = simulation.context.getState(getEnergy=True)
                perturbed_energies[j].append(state.getPotentialEnergy().value_in_unit(kilojoule_per_mole))
        perturbed_energies = np.array(perturbed_energies)
        total_perturbed_energies = np.sum(perturbed_energies, axis=0)
        # compute new total energy by subtracting original value and adding new value
        new_total_energy = np.array(order_parameter_values[:, total_energy_column-1])
        for j, energy_term in enumerate(perturbation):
            new_total_energy -= order_parameter_values[:, energy_term[0]-1]
        new_total_energy += total_perturbed_energies
        all_perturbed_energies.append(new_total_energy)
    all_perturbed_energies = np.array(all_perturbed_energies)
    return all_perturbed_energies.T

def pick_structures(label, conditions_string, metadata_file, reference_structure, num_snapshots=6000, order_parameter_file_name="order_parameters.txt", pdb_trajectory_filename="output.pdb"):
    # parse restrictions (including which file to pull columns from)
    # read in data (including extra files if necessary)
    # go through the data and filter out snapshots that do not satisfy the criteria
    # select a subset of the structures that satisfy the constraints

    def parse_conditions_string(conditions_string):
        conditions = []
        condition_signs = []
        conditions_string = conditions_string.split()
        for condition in conditions_string:
            if "gt" in condition:
                condition_signs.append("+")
                condition = condition.split("gt")
            if "lt" in condition:
                condition_signs.append("-")
                condition = condition.split("lt")
            conditions.append(condition)
        return conditions, condition_signs

    def load_all_data(metadata_file, conditions):
        data_files = list(np.loadtxt(metadata_file, dtype=str)[:,0])
        num_files = len(data_files)
        num_conditions = len(conditions)
        # Load all data into array
        data_array = np.zeros((num_files, num_conditions, num_snapshots))
        for i, data_file in enumerate(data_files):
            all_order_parameters = np.loadtxt(data_file)
            for j, condition in enumerate(conditions):
                data_array[i][j] = all_order_parameters[0:num_snapshots,int(condition[0])-1]

        data_array = np.swapaxes(data_array,1,2)
        return data_array

    def write_selected_pdbs(pdb_files, snapshots_in_pdb_files):
        structure_index = 1
        selected_pdbs = open("%s.pdb" % label, 'w')
        for pdb_file in pdb_files:
            pdb_trajectory_contents = open(pdb_file).read().split("MODEL")[1:]
            pdb_trajectory_contents = ['\n'.join(x.split('\n')[1:]) for x in pdb_trajectory_contents]
            for snapshot in snapshots_in_pdb_files[pdb_files.index(pdb_file)]:
                selected_pdbs.write("MODEL        %d\n" % structure_index)
                selected_pdbs.write(pdb_trajectory_contents[snapshot])
                structure_index += 1
        selected_pdbs.close()

    # Lists
    files_array = list(np.loadtxt(metadata_file, dtype=str)[:,0])
    conditions, condition_signs = parse_conditions_string(conditions_string)
    data_array = load_all_data(metadata_file, conditions)

    # File names and parameters
    output_file_name = label + ".dat"

    structure_index = 1

    # loop over data and output those points that satisfy all conditions
    output_file = open(output_file_name, "w")
    pdb_files = []
    snapshots_in_pdb_files = []

    # loop over files
    for i, data_file in enumerate(files_array):
        # loop over snapshots
        for j, snapshot in enumerate(data_array[i]):
            bad_condition = False
            # loop over conditions
            for k, condition in enumerate(conditions):
                condition_boundary = float(condition[1])
                # If all conditions are satisfied, print out the data
                if condition_signs[k] == "+":
                    if not data_array[i][j][k] > condition_boundary: bad_condition = True
                elif condition_signs[k] == "-":
                    if not data_array[i][j][k] < condition_boundary: bad_condition = True
                else:
                    print("Bad condition argument.")
                    sys.exit()
            if not bad_condition:
                output_file.write("%d %s %s\n" % (structure_index, data_file, j+1))
                pdb_file = data_file.replace(order_parameter_file_name, pdb_trajectory_filename)
                if not pdb_file in pdb_files:
                    pdb_files.append(pdb_file)
                    snapshots_in_pdb_files.append([])
                snapshots_in_pdb_files[pdb_files.index(pdb_file)].append(j)
                structure_index += 1
    output_file.close()
    write_selected_pdbs(pdb_files, snapshots_in_pdb_files)

    pymol_script = open("%s.pml" % label, 'w')
    pymol_script.write("load %s\n" % reference_structure)
    pymol_script.write("load_traj %s.pdb\n" % label)
    object_name = os.path.basename(reference_structure)[:-4]
    pymol_script.write("intra_fit %s\n" % object_name)
    pymol_script.write("smooth\n")
    pymol_script.close()
