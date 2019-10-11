import pandas
import simtk.openmm
import os
import shutil
# imports for accessibility outside
import functionTerms
functionTerms=functionTerms
# Reads pdb file to a table
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
_AWSEMresidues = ['IPR', 'IGL', 'NGP']
xml = f'{__location__}/awsem.xml'

def parsePDB(pdb_file):
    def pdb_line(line):
        return dict(recname=str(line[0:6]).strip(),
                    serial=int(line[6:11]),
                    name=str(line[12:16]).strip(),
                    altLoc=str(line[16:17]),
                    resname=str(line[17:20]).strip(),
                    chainID=str(line[21:22]),
                    resSeq=int(line[22:26]),
                    iCode=str(line[26:27]),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=float(line[54:60]),
                    tempFactor=float(line[60:66]),
                    element=str(line[76:78]),
                    charge=str(line[78:80]))

    with open(pdb_file, 'r') as pdb:
        lines = []
        for line in pdb:
            if len(line) > 6 and line[:6] in ['ATOM  ', 'HETATM']:
                lines += [pdb_line(line)]
    pdb_atoms = pandas.DataFrame(lines)
    pdb_atoms = pdb_atoms[['recname', 'serial', 'name', 'altLoc',
                           'resname', 'chainID', 'resSeq', 'iCode',
                           'x', 'y', 'z', 'occupancy', 'tempFactor',
                           'element', 'charge']]
    return pdb_atoms


def parseConfigTable(config_section):
    """Parses a section of the configuration file as a table"""

    def readData(config_section, a):
        """Filters comments and returns values as a list"""
        temp = config_section.get(a).split('#')[0].split()
        l = []
        for val in temp:
            val = val.strip()
            try:
                x = int(val)
                l += [x]
            except ValueError:
                try:
                    y = float(val)
                    l += [y]
                except ValueError:
                    l += [val]
        return l

    data = []
    for a in config_section:
        if a == 'name':
            columns = readData(config_section, a)
        elif len(a) > 3 and a[:3] == 'row':
            data += [readData(config_section, a)]
        else:
            print(f'Unexpected row {readData(config_section, a)}')
    return pandas.DataFrame(data, columns=columns)


def copy_parameter_files():
    src = f"{__location__}/parameters"
    dest = '.'
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        if os.path.isfile(full_file_name):
            shutil.copy(full_file_name, dest)

class BaseError(Exception):
    pass


class Protein(object):
    def __init__(self, atoms, sequence, k_awsem=1):
        self.atoms = atoms

        #Include real residue name in atoms
        atoms = self.atoms.copy()
        atoms['chain_res'] = atoms['chainID'].astype(str) + '_' + atoms['resSeq'].astype(str)
        sel = atoms[atoms['resname'].isin(_AWSEMresidues)]
        resix = sel['chain_res'].unique()
        assert len(resix) == len(sequence), \
            f'The number of residues {len(resix)} does not agree with the length of the sequence {len(sequence)}'
        atoms.index = atoms['chain_res']
        for r, s in zip(resix, sequence):
            atoms.loc[r, 'real_resname'] = s
        atoms.index = range(len(atoms))
        self.atoms = atoms

        protein_data = atoms[atoms.resname.isin(_AWSEMresidues)].copy()
        # renumber residues
        resix = (protein_data.chainID + '_' + protein_data.resSeq.astype(str))
        res_unique = resix.unique()
        protein_data['resID'] = resix.replace(dict(zip(res_unique, range(len(res_unique)))))
        # renumber atom types
        atom_types_table = {'N': 'n', 'H': 'h', 'CA': 'ca', 'C': 'c', 'O': 'o', 'CB': 'cb'}
        protein_data['atom_list'] = protein_data['name'].replace(atom_types_table)
        protein_data['idx'] = protein_data.index.astype(int)
        self.protein_data = protein_data
        self.atom_lists = protein_data.pivot(index='resID', columns='atom_list', values='idx').fillna(-1).astype(int)
        self.n = self.atom_lists['n'].tolist()
        self.h = self.atom_lists['h'].tolist()
        self.ca = self.atom_lists['ca'].tolist()
        self.c = self.atom_lists['c'].tolist()
        self.o = self.atom_lists['o'].tolist()
        self.cb = self.atom_lists['cb'].tolist()
        self.nres = len(self.atom_lists)
        self.k_awsem = k_awsem
        self.res_type = [r.iloc[0]['resname'] for i, r in protein_data.groupby('resID')]
        self.chain_starts = [c.iloc[0].resID for i, c in protein_data.groupby('chainID')]
        self.chain_ends = [c.iloc[-1].resID for i, c in protein_data.groupby('chainID')]
        self.natoms = len(atoms)
        self.bonds = self._setup_bonds()
        self.seq = sequence
        self.resi = pandas.merge(self.atoms, self.protein_data, how='left').resID.fillna(-1).astype(int).tolist()
        pass

    def _setup_bonds(self):
        bonds = []
        for i in range(self.nres):
            bonds.append((self.ca[i], self.o[i]))
            if not self.res_type[i] == "IGL":
                bonds.append((self.ca[i], self.cb[i]))
            if i not in self.chain_ends:
                bonds.append((self.ca[i], self.ca[i + 1]))
                bonds.append((self.o[i], self.ca[i + 1]))

        for i in range(self.nres):
            if i not in self.chain_starts and not self.res_type[i] == "IGL":
                bonds.append((self.n[i], self.cb[i]))
            if i not in self.chain_ends and not self.res_type[i] == "IGL":
                bonds.append((self.c[i], self.cb[i]))
            if i not in self.chain_starts and i not in self.chain_ends:
                bonds.append((self.n[i], self.c[i]))
        return bonds

    def setup_virtual_sites(self, system, ):
        # set virtual sites
        for i in range(self.nres):
            if i not in self.chain_starts:
                n_virtual_site = simtk.openmm.ThreeParticleAverageSite(self.ca[i - 1], self.ca[i], self.o[i - 1],
                                                                       0.48318, 0.70328, -0.18643)
                system.setVirtualSite(self.n[i], n_virtual_site)
                if not self.res_type[i] == "IPR":
                    h_virtual_site = simtk.openmm.ThreeParticleAverageSite(self.ca[i - 1], self.ca[i], self.o[i - 1],
                                                                           0.84100, 0.89296, -0.73389)
                    system.setVirtualSite(self.h[i], h_virtual_site)
            if i not in self.chain_ends:
                c_virtual_site = simtk.openmm.ThreeParticleAverageSite(self.ca[i], self.ca[i + 1], self.o[i],
                                                                       0.44365, 0.23520, 0.32115)
                # print("Virtual", c[i])
                system.setVirtualSite(self.c[i], c_virtual_site)

    @classmethod
    def fromPDB(cls, pdb, pdbout='CoarseProtein.pdb'):
        """ Initializes a protein form a pdb, making all the atoms coarse-grained"""
        pass

    @classmethod
    def fromCoarsePDB(cls, pdb_file, sequence):
        """ Initializes the protein from an already coarse grained pdb"""
        atoms = parsePDB(pdb_file)
        return cls(atoms, sequence)

    def parseConfigurationFile(self):
        """ Parses the AWSEM configuration file to use for the topology and to set the forces"""
        pass

    def computeTopology(self):
        """ Compute the bonds and angles from the pdb"""
        pass

    @staticmethod
    def CoarseGrain(pdb_table):
        """ Selects AWSEM atoms from a pdb table and returns a table containing only the coarse-grained atoms for AWSEM """
        protein_residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                            'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                            'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                            'SER', 'THR', 'TRP', 'TYR', 'VAL']
        awsem_atoms = ["N", "H", "CA", "C", "O", "CB"]

        # Select coarse grained atoms
        selection = pdb_table[pdb_table.resname.isin(protein_residues) & pdb_table.name.isin(awsem_atoms)].copy()

        # Remove virtual atoms at the end or begining of the chain
        drop_list = []
        for chain in selection.chainID.unique():
            sel = selection[selection.chainID == chain]
            drop_list += list(sel[(sel.resSeq == sel.resSeq.min()) & sel['name'].isin(['N', 'H'])].index)
            drop_list += list(sel[(sel.resSeq == sel.resSeq.max()) & sel['name'].isin(['C'])].index)
        selection = selection.drop(drop_list)

        # Replace resnames
        selection['real_resname'] = selection.resname.copy()
        resname = selection.resname.copy()
        resname[:] = 'NGP'
        resname[selection.resname == 'PRO'] = 'IPR'
        resname[selection.resname == 'GLY'] = 'IGL'
        selection.resname = resname

        # CB element is B
        selection.loc[selection['name'] == 'CB', 'element'] = 'B'

        # Reorder atoms
        selection.name = pandas.Categorical(selection.name, awsem_atoms)
        selection.sort_values(['chainID', 'resSeq', 'name'])

        # Prepare virtual sites
        for c, chain in selection.groupby('chainID'):
            first = chain.resSeq.min()
            last = chain.resSeq.max()
            for i, residue in chain.groupby('resSeq'):
                idx = dict(zip(residue.name, residue.index))
                pos = dict(zip(residue.name, [residue.loc[i, ['x', 'y', 'z']] for i in residue.index]))

                if i != first:
                    if 'N' in idx.keys():
                        selection.loc[idx['N'], ['x', 'y', 'z']] = 0.48318 * pos_im['CA'] + 0.70328 * pos[
                            'CA'] - 0.18643 * pos_im['O']
                    if 'H' in idx.keys():
                        selection.loc[idx['H'], ['x', 'y', 'z']] = 0.84100 * pos_im['CA'] + 0.89296 * pos[
                            'CA'] - 0.73389 * pos_im['O']
                    if 'C' in idx.keys():
                        selection.loc[idx_im['C'], ['x', 'y', 'z']] = 0.44365 * pos_im['CA'] + 0.23520 * pos[
                            'CA'] + 0.32115 * pos_im['O']

                pos_im = pos.copy()
                idx_im = idx.copy()
        # Renumber
        selection['serial'] = range(len(selection))
        return selection

    @staticmethod
    def write_sequence(Coarse, seq_file='protein.seq'):
        from Bio.PDB.Polypeptide import three_to_one
        protein_data = Coarse[Coarse.resname.isin(_AWSEMresidues)].copy()
        resix = (protein_data.chainID + '_' + protein_data.resSeq.astype(str))
        res_unique = resix.unique()
        protein_data['resID'] = resix.replace(dict(zip(res_unique, range(len(res_unique)))))
        protein_sequence = [r.iloc[0]['real_resname'] for i, r in protein_data.groupby('resID')]
        protein_sequence_one = [three_to_one(a) for a in protein_sequence]

        with open(seq_file, 'w+') as ps:
            ps.write(''.join(protein_sequence_one))

def addNonBondedExclusions(oa, force):
    cb_fixed = [x if x > 0 else y for x, y in zip(oa.cb, oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    for e1 in none_cb_fixed:
        for e2 in none_cb_fixed:
            if e1 > e2:
                continue
            force.addExclusion(e1, e2)
    for e1 in none_cb_fixed:
        for e2 in cb_fixed:
            force.addExclusion(e1, e2)
    for e1 in none_cb_fixed:
        for e2 in none_cb_fixed:
            if e1 > e2:
                continue
            force.addExclusion(e1, e2)
    for e1 in none_cb_fixed:
        for e2 in cb_fixed:
            force.addExclusion(e1, e2)


def Forces():
    pass


def test_Protein_fromCoarsePDB():
    pass


def test_Protein_fromPDB():
    pass


def test_Protein_parseConfigurationFile():
    pass


def test_Protein_computeTopology():
    pass
