from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from pdbfixer import *
import mdtraj as md
from Bio.PDB.Polypeptide import *
from Bio.PDB.PDBParser import PDBParser
from itertools import combinations
import numpy as np

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
		self.resi = [x.residue.index for x in list(self.pdb.topology.atoms())]
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
	k_con = 40
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
	k_chain = 20
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
	k_chi = 80
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

	excl.setCutoffDistance(r_excl)

	excl.createExclusionsFromBonds(bonds, 1)
	system.addForce(excl)

def apply_rama_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type
	# add Rama potential
	# Still need to add proline parameters and secondary structure biases
	k_rama = 50
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

def apply_direct_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds
	# add direct contact
	# Still need to add residue specific parameters
	k_direct = 0.00
	gamma = 1
	r_min = .45
	r_max = .65
	eta = 10
	direct = CustomNonbondedForce("-k_direct*gamma*theta; theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r))); gamma=%f; eta=%f" % (gamma, eta))
	direct.addGlobalParameter("k_direct", k_direct)
	direct.addGlobalParameter("rmin", r_min)
	direct.addGlobalParameter("rmax", r_max)
	for i in range(natoms):
		direct.addParticle()
	direct.addInteractionGroup([x for x in cb if x > 0], [x for x in cb if x > 0])

	direct.createExclusionsFromBonds(bonds, 9)
	system.addForce(direct)

def apply_mediated_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi
	# add mediated contact
	# Still need to add residue specific parameters
	k_mediated = 0.001
	gamma_water = 0
	gamma_protein = 1
	r_min = .65
	r_max = .95
	eta = 10
	density_r_min = 0.45
	density_r_max = 0.65
	density_threshold = 2.6
	interaction_cutoff = .1
	min_sequence_separation = 10
	mediated = CustomGBForce()
	mediated.addGlobalParameter("k_mediated", k_mediated)
	mediated.addGlobalParameter("gamma_water", gamma_water)
	mediated.addGlobalParameter("gamma_protein", gamma_protein)
	mediated.addPerParticleParameter("index")
	include_pairwise_interaction = [0]*natoms*natoms
	for i in range(natoms):
		for j in range(natoms):
			# if this is not a cbeta atom, don't include in pairwise interaction
			if i in cb and j in cb and abs(resi[i]-resi[j]) >= min_sequence_separation:
				include_pairwise_interaction[i+j*natoms] = 1

	mediated.addTabulatedFunction("include_pairwise_interaction", Discrete2DFunction(natoms, natoms, include_pairwise_interaction))
	mediated.addComputedValue("rho", "0.25*(1+tanh(eta*(r-%f)))*(1+tanh(eta*(%f-r))); eta=%f" % (density_r_min, density_r_max, eta), CustomGBForce.ParticlePair)
	sigma_water="0.25*(1-tanh(%f*(rho1-%f)))*(1-tanh(%f*(rho2-%f)))" % (eta, density_threshold, eta, density_threshold)
	mediated.addEnergyTerm("-k_mediated*include_pairwise_interaction(index1,index2)*theta*(gamma_water*%s+gamma_protein*(1-%s)); theta=0.25*(1+tanh(eta*(r-%f)))*(1+tanh(eta*(%f-r))); gamma_water=%f; gamma_protein=%f; eta=%f" % (sigma_water, sigma_water, r_min, r_max, gamma_water, gamma_protein, eta), CustomGBForce.ParticlePair)

	# set interaction cutoff
	mediated.setCutoffDistance(interaction_cutoff)

	# add all particles to the force
	for i in range(natoms):
		mediated.addParticle([i])

	# find all pairs to include in the force
	interactions = []
	for cbi, cbj in combinations([x for x in cb if x >= 0], 2):
		if abs(resi[cbi]-resi[cbj]) >= 1:
			interactions.append((cbi, cbj))

	# get list of exclusions
	exclusions = get_exclusions(oa, interactions)

	# apply exclusions
	for exclusion in exclusions:
		mediated.addExclusion(exclusion[0], exclusion[1])
	
	system.addForce(mediated)

def get_exclusions(oa, interactions):
	natoms = oa.natoms
	all_interactions = combinations(range(natoms), 2)
	exclusions = [x for x in all_interactions if x not in interactions]
	return exclusions

def apply_contact_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi

	# setup parameters
	k_direct = 0.01
	k_mediated = 0.01
	k_burial = 1

	r_direct_min = .45
	r_direct_max = .65

	r_mediated_min = .65
	r_mediated_max = .95

	r_density_min = 0.45
	r_density_max = 0.65
	rho0 = 2.6

	interaction_cutoff_distance = .1
	min_sequence_separation = 10

	eta_direct = 50
	eta_mediated = 50
	eta_sigma = 70
	eta_burial = 40

	rho_min_1 = 0
	rho_min_2 = 3
	rho_min_3 = 6
	rho_max_1 = 3
	rho_max_2 = 6
	rho_max_3 = 9

	# setup matrix of pairwise weights of interactions
	# this also serves to exclude interactions below the
	# minimum sequence separation by setting those gamma values to zero
	gamma_direct = [0]*natoms*natoms
	gamma_water = [0]*natoms*natoms
	gamma_protein = [0]*natoms*natoms
	gamma_burial_1 = [0]*natoms
	gamma_burial_2 = [0]*natoms
	gamma_burial_3 = [0]*natoms
	for i in range(natoms):
		if i in cb:
			gamma_burial_1[i] = 1
			gamma_burial_2[i] = 1
			gamma_burial_3[i] = 1
		for j in range(natoms):
			# if this is not a cbeta atom, don't include in pairwise interaction
			if i in cb and j in cb and abs(resi[i]-resi[j]) >= min_sequence_separation:
				gamma_direct[i+j*natoms] = 1
				gamma_water[i+j*natoms] = 1
				gamma_protein[i+j*natoms] = 1

	# create contact force
	contact = CustomGBForce()
	# add per-particle parameters
	contact.addPerParticleParameter("index")
	# add global parameters
	contact.addGlobalParameter("k_direct", k_direct)
	contact.addGlobalParameter("eta_direct", eta_direct)
	contact.addGlobalParameter("k_mediated", k_mediated)
	contact.addGlobalParameter("eta_mediated", eta_mediated)
	contact.addGlobalParameter("k_burial", k_burial)
	contact.addGlobalParameter("rho0", rho0)
	contact.addGlobalParameter("eta_burial", eta_burial)
	contact.addGlobalParameter("eta_sigma", eta_sigma)
	contact.addGlobalParameter("r_direct_min", r_direct_min)
	contact.addGlobalParameter("r_direct_max", r_direct_max)
	contact.addGlobalParameter("r_mediated_min", r_mediated_min)
	contact.addGlobalParameter("r_mediated_max", r_mediated_max)
	contact.addGlobalParameter("r_density_min", r_density_min)
	contact.addGlobalParameter("r_density_max", r_density_max)
	contact.addGlobalParameter("rho_min_1", rho_min_1)
	contact.addGlobalParameter("rho_min_2", rho_min_2)
	contact.addGlobalParameter("rho_min_3", rho_min_3)
	contact.addGlobalParameter("rho_max_1", rho_max_1)
	contact.addGlobalParameter("rho_max_2", rho_max_2)
	contact.addGlobalParameter("rho_max_3", rho_max_3)
	# add tabulated functions
	contact.addTabulatedFunction("gamma_direct", Discrete2DFunction(natoms, natoms, gamma_direct))
	contact.addTabulatedFunction("gamma_water", Discrete2DFunction(natoms, natoms, gamma_water))
	contact.addTabulatedFunction("gamma_protein", Discrete2DFunction(natoms, natoms, gamma_protein))
	contact.addTabulatedFunction("gamma_burial_1", Discrete1DFunction(gamma_burial_1))
	contact.addTabulatedFunction("gamma_burial_2", Discrete1DFunction(gamma_burial_2))
	contact.addTabulatedFunction("gamma_burial_3", Discrete1DFunction(gamma_burial_3))

	# compute the density, exclusions will be taken care of
	# using the addExclusion function below
	contact.addComputedValue("rho", "0.25*(1+tanh(eta_direct*(r-r_density_min)))*(1+tanh(eta_direct*(r_density_max-r)))", CustomGBForce.ParticlePair)
	# add direct interaction
	contact.addEnergyTerm("-k_direct*gamma_direct(index1,index2)*theta_direct; theta_direct=0.25*(1+tanh(eta_direct*(r-r_direct_min)))*(1+tanh(eta_direct*(r_direct_max-r)))", CustomGBForce.ParticlePair)
	# add mediated interaction
	contact.addEnergyTerm("-k_mediated*theta_mediated*(gamma_water(index1,index2)*sigma_water+gamma_protein(index1,index2)*(1-sigma_water)); theta_mediated=0.25*(1+tanh(eta_mediated*(r-r_mediated_min)))*(1+tanh(eta_mediated*(r_mediated_max-r))); sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho0)))*(1-tanh(eta_sigma*(rho2-rho0)))", CustomGBForce.ParticlePair)
	# add burial term
	burial_string = ''.join(["-0.5*k_burial*gamma_burial_%d(index)*(tanh(eta_burial*(rho-rho_min_%d))+tanh(eta_burial*(rho_max_%d-rho)))+" % (x, x, x) for x in range(1,4)])[:-1]
	contact.addEnergyTerm(burial_string, CustomGBForce.SingleParticle)

	# set interaction cutoff
	contact.setCutoffDistance(interaction_cutoff_distance)

	# add all particles to the force
	for i in range(natoms):
		contact.addParticle([i])

	# find all pairs to include in the force
	interactions = []
	for cbi, cbj in combinations([x for x in cb if x >= 0], 2):
		if abs(resi[cbi]-resi[cbj]) >= 1:
			interactions.append((cbi, cbj))

	# get list of exclusions
	exclusions = get_exclusions(oa, interactions)

	# apply exclusions
	for exclusion in exclusions:
		contact.addExclusion(exclusion[0], exclusion[1])
	
	system.addForce(contact)

def apply_beta_term(oa):
	system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi

	def lambda_1(i, j):
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 1.37
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 1.36
		elif abs(j-i) >= 45:
			return 1.17
	def lambda_2(i, j):
		extra_terms = -0.5*alpha_1(i,j)*np.log(p_hb(i, j))-0.25*alpha_2(i, j)*(np.log(p_nhb(i+1,j-1)+np.log(p_nhb(i-1,j+1))))-alpha_3(i, j)*(np.log(p_anti(i))+np.log(p_anti(j)))
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 3.89+extra_terms
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 3.50+extra_terms
		elif abs(j-i) >= 45:
			return 3.52+extra_terms
	def lambda_3(i, j):
		extra_terms = -alpha_4(i,j)*np.log(p_parhb(i+1,j))-alpha_5(i,j)*np.log(p_par(i+1))+alpha_4*np.log(p_par(j))
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 0.0+extra_terms
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 3.47+extra_terms
		elif abs(j-i) >= 45:
			return 3.62+extra_terms
	def alpha_1(i,j):
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 1.3
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 1.3
		elif abs(j-i) >= 45:
			return 1.3
	def alpha_2(i,j):
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 1.32
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 1.32
		elif abs(j-i) >= 45:
			return 1.32
	def alpha_3(i,j):
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 1.22
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 1.22
		elif abs(j-i) >= 45:
			return 1.22	
	def alpha_4(i,j):
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 0.0
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 0.33
		elif abs(j-i) >= 45:
			return 0.33
	def alpha_5(i,j):
		if abs(j-i) >= 4 and abs(j-i) < 18:
			return 0.0
		elif abs(j-i) >= 18 and abs(j-i) < 45:
			return 1.01
		elif abs(j-i) >= 45:
			return 1.01
	
	# add beta potential
	# setup parameters
	k_beta = 1.0
	lambda_1 = [0]*nres*nres
	for i in range(nres):
		for j in range(nres):
			lambda_1[i+j*nres] = lambda_1(i,j)
			lambda_2[i+j*nres] = lambda_2(i,j)
			lambda_3[i+j*nres] = lambda_3(i,j)

	r_ON = .298
	sigma_NO = .068
	r_OH = .206
	sigma_HO = .076
	eta_beta_1 = 10.0
	eta_beta_2 = 5.0
	r_HB_c = 1.2

	theta_ij =   "exp(-(r_Oi_Nj-r_ON)^2/(2*sigma_NO^2)-(r_Oi_Hj-r_OH)^2/(2*sigma_HO^2))"
	theta_ji =   "exp(-(r_Oj_Ni-r_ON)^2/(2*sigma_NO^2)-(r_Oj_Hi-r_OH)^2/(2*sigma_HO^2))"
	theta_jip2 = "exp(-(r_Oj_Nip2-r_ON)^2/(2*sigma_NO^2)-(r_Oj_Hip2-r_OH)^2/(2*sigma_HO^2))"
	nu_i = "0.5*(1+tanh(eta_beta_1*(r_CAim2_CAip2-r_HB_c)))"
	nu_j = "0.5*(1+tanh(eta_beta_2*(r_CAjm2_CAjp2-r_HB_c)))"

	# Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
	# 1  2  3  4     5     6     7
	beta_string_1 = "-k_beta*lambda_1(index_i,index_j)*theta_ij*nu_i*nu_j;theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
					nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)" % (theta_ij, nu_i, nu_j)

	# Oi Nj Hj Oj Ni Hi CAi-2 CAi+2 CAj-2 CAj+2
	# 1  2  3  4  5  6  7     8     9     10
	beta_string_2 = "-k_beta*lambda_2(index_i,index_j)*theta_ij*theta_ji*nu_i*nu_j;\
					theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
					theta_ji=%s;r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
					nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)" % (theta_ij, theta_ji, nu_i, nu_j)

	# Oi Nj Hj Oj Ni+2 Hi+2 CAi-2 CAi+2 CAj-2 CAj+2
	# 1  2  3  4  5    6    7     8     9     10
	beta_string_3 = "-k_beta*lambda_3(index_i,index_j)*theta_ij*theta_jip2*nu_i*nu_j;\
					theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
					theta_ji=%s;r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
					theta_jip2=%s;r_Oj_Nip2=distance(p4,p5);r_Oj_Hip2=distance(p4,p6);\
					nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)" % (theta_ij, theta_ji, theta_jip2, nu_i, nu_j)

	beta_1 = CustomCompoundBondForce(7, beta_string_1)
	beta_2 = CustomCompoundBondForce(10, beta_string_2)
	beta_3 = CustomCompoundBondForce(10, beta_string_3)
	# add parameters to force
	for i in range(1,4):
		locals()["beta_%d" % i].addGlobalParameter("k_beta", k_beta)
		locals()["beta_%d" % i].addGlobalParameter("r_ON", r_ON)
		locals()["beta_%d" % i].addGlobalParameter("sigma_NO", sigma_NO)
		locals()["beta_%d" % i].addGlobalParameter("r_OH", r_OH)
		locals()["beta_%d" % i].addGlobalParameter("sigma_HO", sigma_HO)
		locals()["beta_%d" % i].addGlobalParameter("eta_beta_1", eta_beta_1)
		locals()["beta_%d" % i].addGlobalParameter("eta_beta_2", eta_beta_2)
		locals()["beta_%d" % i].addGlobalParameter("r_HB_c", r_HB_c)
		locals()["beta_%d" % i].addPerBondParameter("index_i")
		locals()["beta_%d" % i].addPerBondParameter("index_j")
	beta_1.addTabulatedFunction("lambda_1", Discrete2DFunction(nres, nres, lambda_1))
	beta_2.addTabulatedFunction("lambda_2", Discrete2DFunction(nres, nres, lambda_2))
	beta_3.addTabulatedFunction("lambda_3", Discrete2DFunction(nres, nres, lambda_3))

	for i in range(nres):
		for j in range(i, nres):
			if i-2 < 0 or i+2 >= nres or \
			   j-2 < 0 or j+2 >= nres:
			   continue
			if not res_type[j] == "IPR":
				beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
			if not res_type[i] == "IPR" and not res_type[j] == "IPR":
				beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
			if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
				beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])

	system.addForce(beta_1)
	system.addForce(beta_2)
	system.addForce(beta_3)