from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from pdbfixer import *
import mdtraj as md
from Bio.PDB.Polypeptide import *
from Bio.PDB.PDBParser import PDBParser

def identify_terminal_residues(pdb_filename):
	# identify terminal residues
	p = PDBParser()
	structure = p.get_structure('X', pdb_filename)
	terminal_residues = {}
	for model in structure:
		for chain in model:
			residues = list(chain.get_residues())
			terminal_residues[chain.id] = (residues[0].id[1], residues[-1].id[1])

	return terminal_residues

def prepare_pdb(pdb_filename, chains_to_simulate, input_pdb_filename=None):
	# for more information about PDB Fixer, see:
	# http://htmlpreview.github.io/?https://raw.github.com/pandegroup/pdbfixer/master/Manual.html
	# fix up input pdb
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
	PDBFile.writeFile(fixer.topology, fixer.positions, open('pdbfixeroutput.pdb', 'w'))

	# identify terminal residues
	terminal_residues = identify_terminal_residues('pdbfixeroutput.pdb')

	# process pdb for input into OpenMM
	if input_pdb_filename == None:
		input_pdb_filename = pdb_filename.split('.')[0] + '-openmmawsem.pdb'
	output = open(input_pdb_filename, 'w')
	for line in open("pdbfixeroutput.pdb"):
		splitline = line.split()
		if splitline[0] == "ATOM":
			_, atom_index, atom_type, res_type, chain, res_index, x, y, z, _, _, element = splitline
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
			if splitline[3] == "GLY":
				line[17:20] = "IGL"
			elif splitline[3] == "PRO":
				line[17:20] = "IPR"
			else:
				line[17:20] = "NGP"
			if atom_type == "CB":
				line[77] = "B"
			line=''.join(line)    
			output.write(line)
	output.close()

def build_lists_of_atoms(nres, residues):
	# build lists of atoms, residue types, and bonds
	n, h, ca, c, o, cb, res_type = [], [], [], [], [], [], []
	for residue in residues:
		res_type.append(residue.name)
		atom_lists = [n, h, ca, c, o, cb]
		residue_atoms = [x.index for x in residue.atoms()]
		if residue.index == 0:
			atom_lists.remove(n)
			n.append(-1)
			atom_lists.remove(h)
			h.append(-1)
		if residue.index == nres-1:
			atom_lists.remove(c)
			c.append(-1)
		if residue.name == "IGL" and cb in atom_lists:
			atom_lists.remove(cb)
			cb.append(-1)
		if residue.name == "IPR" and h in atom_lists:
			atom_lists.remove(h)
			h.append(-1)
		for atom, atom_list in zip(residue_atoms, atom_lists):
				atom_list.append(atom)
	return n, h, ca, c, o, cb, res_type

def setup_virtual_sites(nres, system, n, h, ca, c, o, cb, res_type):
	# set virtual sites
	for i in range(nres):
		if i > 0:
			n_virtual_site = ThreeParticleAverageSite(ca[i-1], ca[i], o[i-1], 0.48318, 0.70328, -0.18643)
			system.setVirtualSite(n[i], n_virtual_site)
			if not res_type[i] == "IPR":
				h_virtual_site = ThreeParticleAverageSite(ca[i-1], ca[i], o[i-1], 0.84100, 0.89296, -0.73389)
				system.setVirtualSite(h[i], h_virtual_site)
		if  i+1 < nres:
			c_virtual_site = ThreeParticleAverageSite(ca[i], ca[i+1], o[i], 0.44365, 0.23520, 0.32115)
			system.setVirtualSite(c[i], c_virtual_site)

def setup_bonds(nres, n, h, ca, c, o, cb, res_type):
	bonds = []
	for i in range(nres):
		bonds.append((ca[i], o[i]))
		if not res_type[i] == "IGL":
			bonds.append((ca[i], cb[i]))
		if  i+1 < nres:            
			bonds.append((ca[i], ca[i+1]))
			bonds.append((o[i], ca[i+1]))

	for i in range(nres):      
		if not i == 0 and not res_type[i] == "IGL":
			bonds.append((n[i], cb[i]))
		if not i+1 == nres and not res_type[i] == "IGL":
			bonds.append((c[i], cb[i]))
		if not i == 0 and not i+1 == nres:
			bonds.append((n[i], c[i]))
	return bonds

class OpenMMAWSEMSystem:
	def __init__(self, pdb_filename, xml_filename='awsem.xml'):
		# read PDB
		self.pdb = PDBFile(pdb_filename)
		self.forcefield = ForceField(xml_filename)
		self.system = self.forcefield.createSystem(self.pdb.topology)
		# define convenience variables
		self.nres = self.pdb.topology.getNumResidues()
		self.natoms = self.pdb.topology.getNumAtoms()
		self.residues = list(self.pdb.topology.residues())
		# build lists of atoms and residue types
		self.n, self.h, self.ca, self.c, self.o, self.cb, self.res_type = build_lists_of_atoms(self.nres, self.residues)
		# setup virtual sites
		setup_virtual_sites(self.nres, self.system, self.n, self.h, self.ca, self.c, self.o, self.cb, self.res_type)
		# setup bonds
		self.bonds = setup_bonds(self.nres, self.n, self.h, self.ca, self.c, self.o, self.cb, self.res_type)
		# identify terminal_residues
		self.terminal_residues = identify_terminal_residues(pdb_filename)

def apply_con_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type
	# add con forces
	con = HarmonicBondForce()
	k_con = 400
	for i in range(nres):
		con.addBond(ca[i], o[i], .243, k_con)
		if not res_type[i] == "IGL":
			con.addBond(ca[i], cb[i], .154, k_con)
		if  i+1 < nres:            
			con.addBond(ca[i], ca[i+1], .380, k_con)
			con.addBond(o[i], ca[i+1], .282, k_con)

	system.addForce(con)

def apply_chain_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type
	# add chain forces
	chain = HarmonicBondForce()
	k_chain = 200
	for i in range(nres):      
		if not i == 0 and not res_type[i] == "IGL":
			chain.addBond(n[i], cb[i], .246, k_chain)
		if not i+1 == nres and not res_type[i] == "IGL":
			chain.addBond(c[i], cb[i], .270, k_chain)
		if not i == 0 and not i+1 == nres:
			chain.addBond(n[i], c[i],  .246, k_chain)

	system.addForce(chain)

def apply_chi_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type
	# add chi forces 
	# The sign of the equilibrium value is opposite and magnitude differs slightly
	k_chi = 800
	chi0 = .0093
	chi = CustomCompoundBondForce(4, "k_chi*(chi-chi0)^2;\
								  chi=crossproduct_x*r_cacb_x+crossproduct_y*r_cacb_y+crossproduct_z*r_cacb_z;\
								  crossproduct_x=(u2*v3-u3*v2);\
								  crossproduct_y=(u3*v1-u1*v3);\
								  crossproduct_z=(u1*v2-u2*v1);\
								  r_cacb_x=x4-x1;\
								  r_cacb_y=y4-y1;\
								  r_cacb_z=z4-z1;\
								  u1=x1-x2; u2=y1-y2; u3=z1-z2;\
								  v1=x3-x2; v2=y3-y2; v3=z3-z2")
	chi.addGlobalParameter("k_chi", k_chi)
	chi.addGlobalParameter("chi0", chi0)
	for i in range(nres):
		if not i == 0 and not i+1 == nres and not res_type[i] == "IGL":
			chi.addBond([ca[i], c[i], n[i], cb[i]])

	system.addForce(chi)

def apply_excl_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds
	# add excluded volume
	# Still need to add element specific parameters
	k_excl = 80
	r_excl = .35
	excl = CustomNonbondedForce("k_excl*step(r0-r)*(r-r0)^2")
	excl.addGlobalParameter("k_excl", k_excl)
	excl.addGlobalParameter("r0", r_excl)
	for i in range(natoms):
		excl.addParticle()
	excl.addInteractionGroup(ca, ca)
	excl.addInteractionGroup([x for x in cb if x > 0], [x for x in cb if x > 0])
	excl.addInteractionGroup(o, o)

	excl.createExclusionsFromBonds(bonds, 1)
	system.addForce(excl)

def apply_rama_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type
	# add Rama potential
	# Still need to add proline parameters and secondary structure biases
	k_rama = 5
	num_rama_wells = 3
	w = [1.3149, 1.32016, 1.0264]
	sigma = [15.398, 49.0521, 49.0954]
	omega_phi = [0.15, 0.25, 0.65]
	phi_i = [-1.74, -1.265, 1.041]
	omega_psi = [0.65, 0.45, 0.25]
	psi_i = [2.138, 0.318, 0.78]
	rama_function = ''.join(["w%d*exp(-sigma%d*(omega_phi%d*phi_term%d^2+omega_psi%d*psi_term%d^2))+" \
							  % (i, i, i, i, i, i) for i in range(num_rama_wells)])[:-1]+';'
	rama_parameters = ''.join(["w%d=%f; sigma%d=%f; omega_phi%d = %f;\
							   phi_term%d=cos(phi_i-%f)-1; phi_i=dihedral(p1, p2, p3, p4);\
							   omega_psi%d=%f;\
							   psi_term%d=cos(psi_i-%f)-1; psi_i=dihedral(p2, p3, p4, p5);" \
							  % (i, w[i], i, sigma[i], i, omega_phi[i], i, phi_i[i], i, omega_psi[i],\
								i, psi_i[i]) for i in range(num_rama_wells)])[:-1]
	rama_string = rama_function+rama_parameters
	rama = CustomCompoundBondForce(5, rama_string)

	for i in range(nres):
		if not i == 0 and not i+1 == nres and not res_type[i] == "IGL":
			rama.addBond([c[i-1], n[i], ca[i], c[i], n[i+1]])

	system.addForce(rama)

