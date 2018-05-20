from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from pdbfixer import *
import mdtraj as md
from Bio.PDB.Polypeptide import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from itertools import combinations
import numpy as np

se_map_3_letter = dict(zip(("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"), (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)))

se_map_1_letter = dict(zip(("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"), (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)))

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

    #Read sequence
    structure = PDBParser().get_structure('pdbfixeroutput','pdbfixeroutput.pdb')
    res_names = list([x.resname for x in structure.get_residues()])

    # identify terminal residues
    terminal_residues = identify_terminal_residues('pdbfixeroutput.pdb')

    # process pdb for input into OpenMM
    #Selects only atoms needed for the awsem topology
    if input_pdb_filename == None:
        input_pdb_filename = pdb_filename.split('.')[0] + '-openmmawsem.pdb'

    output = open(input_pdb_filename, 'w')
    counter=0
    for line in open("pdbfixeroutput.pdb"):
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

    return res_names

def prepare_virtual_sites(pdb_file):
    p = PDBParser(QUIET=True)
    structure=p.get_structure('X',pdb_file,)
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
        #print(residue_atoms)
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
                #print(atype,atom)
                atom_lists[atype].append(atom)
    #[print(key,len(atom_lists[key])) for key in atom_lists]

    return atom_lists, res_types

def setup_virtual_sites(nres, system, n, h, ca, c, o, cb, res_type):
    # set virtual sites
    for i in range(nres):
        if i > 0:
            n_virtual_site = ThreeParticleAverageSite(ca[i-1], ca[i], o[i-1],
                                                      0.48318, 0.70328, -0.18643)
            system.setVirtualSite(n[i], n_virtual_site)
            if not res_type[i] == "IPR":
                h_virtual_site = ThreeParticleAverageSite(ca[i-1], ca[i], o[i-1],
                                                          0.84100, 0.89296, -0.73389)
                system.setVirtualSite(h[i], h_virtual_site)
        if  i+1 < nres:
            c_virtual_site = ThreeParticleAverageSite(ca[i], ca[i+1], o[i],
                                                      0.44365, 0.23520, 0.32115)
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
    def __init__(self, pdb_filename, res_names, xml_filename='awsem.xml'):
        # get full residue names
        self.res_names = res_names
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
        self.atom_lists,self.res_type=build_lists_of_atoms(self.nres, self.residues)
        #print(str(self.atom_lists))
        self.n =self.atom_lists['n']
        self.h =self.atom_lists['h']
        self.ca=self.atom_lists['ca']
        self.c =self.atom_lists['c']
        self.o =self.atom_lists['o']
        self.cb=self.atom_lists['cb']
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
    k_con = 60
    for i in range(nres):
        con.addBond(ca[i], o[i], .243, k_con)
        if not res_type[i] == "IGL":
            con.addBond(ca[i], cb[i], .154, k_con)
        if  i+1 < nres:
            con.addBond(ca[i], ca[i+1], .380, k_con)
            con.addBond(o[i], ca[i+1], .282, k_con)

    system.addForce(con)

def apply_chain_term(oa):
    # add chain forces
    chain = HarmonicBondForce()
    k = np.array([60., 60., 60.])* 2 *4.184 * 100.      # kcal/A^2 to kJ/nm^2
    x = np.array([2.459108, 2.519591, 2.466597])/10. # nm to A
    #x = np.array([2.46, 2.7, 2.46])/10. # nm to A
    #x = np.array([2.46, 2.52, 2.42])/10. # nm to A
    #x = np.array([2.4545970985006895, 2.564555486626491, 2.548508839171672])/10.
    #x = np.array([2.455, 2.565, 2.548])/10.
    for i in range(oa.nres):
        if not i == 0 and not oa.res_type[i] == "IGL":
            chain.addBond(oa.n[i], oa.cb[i], x[0], k[0])
        if not i+1 == oa.nres and not oa.res_type[i] == "IGL":
            chain.addBond(oa.c[i], oa.cb[i], x[1], k[1])
        if not i == 0 and not i+1 == oa.nres:
            chain.addBond(oa.n[i], oa.c[i],  x[2], k[2])
    oa.system.addForce(chain)

def apply_chi_term(oa):
    system, nres, n, h, ca, c, o, cb, res_type = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type
    # add chi forces
    # The sign of the equilibrium value is opposite and magnitude differs slightly
    # k_chi = 80
    # chi0 = .0093
    k_chi = 60 * 4.184 # kCal to kJ
    chi0 = -0.71 # none dimensional
    # chi = CustomCompoundBondForce(4, "k_chi*(chi-chi0)^2;\
    #                               chi=crossproduct_x*r_cacb_x+crossproduct_y*r_cacb_y+crossproduct_z*r_cacb_z;\
    #                               crossproduct_x=(u2*v3-u3*v2);\
    #                               crossproduct_y=(u3*v1-u1*v3);\
    #                               crossproduct_z=(u1*v2-u2*v1);\
    #                               r_cacb_x=x4-x1;\
    #                               r_cacb_y=y4-y1;\
    #                               r_cacb_z=z4-z1;\
    #                               u1=x1-x2; u2=y1-y2; u3=z1-z2;\
    #                               v1=x3-x2; v2=y3-y2; v3=z3-z2")
    chi = CustomCompoundBondForce(4, "k_chi*(chi*norm-chi0)^2;\
                                  chi=crossproduct_x*r_cacb_x+crossproduct_y*r_cacb_y+crossproduct_z*r_cacb_z;\
                                  crossproduct_x=(u_y*v_z-u_z*v_y);\
                                  crossproduct_y=(u_z*v_x-u_x*v_z);\
                                  crossproduct_z=(u_x*v_y-u_y*v_x);\
                                  norm=1/((u_x*u_x+u_y*u_y+u_z*u_z)*(v_x*v_x+v_y*v_y+v_z*v_z)*(r_cacb_x*r_cacb_x+r_cacb_y*r_cacb_y+r_cacb_z*r_cacb_z))^0.5;\
                                  r_cacb_x=x1-x4;\
                                  r_cacb_y=y1-y4;\
                                  r_cacb_z=z1-z4;\
                                  u_x=x1-x2; u_y=y1-y2; u_z=z1-z2;\
                                  v_x=x3-x1; v_y=y3-y1; v_z=z3-z1;")
    # chi = CustomCompoundBondForce(4, "x4")
    # chi = CustomCompoundBondForce(4, "z1")
    chi.addGlobalParameter("k_chi", k_chi)
    chi.addGlobalParameter("chi0", chi0)
    for i in range(nres):
        if not i == 0 and not i+1 == nres and not res_type[i] == "IGL":
            chi.addBond([ca[i], c[i], n[i], cb[i]])
    # i = 1
    # chi.addBond([ca[i], c[i], n[i], cb[i]])
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
    k_rama = -2 * 4.184 # in kJ/mol
    num_rama_wells = 3 # probably better to apply the ssweight bias in another function
    w = [1.3149, 1.32016, 1.0264]
    sigma = [15.398, 49.0521, 49.0954]
    omega_phi = [0.15, 0.25, 0.65]
    phi_i = [-1.74, -1.265, 1.041]
    omega_psi = [0.65, 0.45, 0.25]
    psi_i = [2.138, -0.318, 0.78]
    rama_function = ''.join(["w%d*exp(-sigma%d*(omega_phi%d*phi_term%d^2+omega_psi%d*psi_term%d^2))+" \
                              % (i, i, i, i, i, i) for i in range(num_rama_wells)])[:-1]
    rama_function = 'k_rama*(' + rama_function + ");"
    # rama_parameters = ''.join(["w%d=%f; sigma%d=%f; omega_phi%d = %f;\
    #                            phi_term%d=cos(phi_i-%f)-1; phi_i=dihedral(p1, p2, p3, p4);\
    #                            omega_psi%d=%f;\
    #                            psi_term%d=cos(psi_i-%f)-1; psi_i=dihedral(p2, p3, p4, p5);" \
    #                           % (i, w[i], i, sigma[i], i, omega_phi[i], i, phi_i[i], i, omega_psi[i],\
    #                             i, psi_i[i]) for i in range(num_rama_wells)])[:-1]
    rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                            psi_term{i}=cos(psi_{i}-psi0{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                             for i in range(num_rama_wells)])
    rama_string = rama_function+rama_parameters

    # rama_string = "dihedral(p1, p2, p3, p4);"
    rama = CustomCompoundBondForce(5, rama_string)
    for i in range(num_rama_wells):
        rama.addGlobalParameter(f"k_rama", k_rama)
        rama.addGlobalParameter(f"w{i}", w[i])
        rama.addGlobalParameter(f"sigma{i}", sigma[i])
        rama.addGlobalParameter(f"omega_phi{i}", omega_phi[i])
        rama.addGlobalParameter(f"omega_psi{i}", omega_psi[i])
        rama.addGlobalParameter(f"phi0{i}", phi_i[i])
        rama.addGlobalParameter(f"psi0{i}", psi_i[i])
    for i in range(nres):
        if not i == 0 and not i+1 == nres and not res_type[i] == "IGL":
            rama.addBond([c[i-1], n[i], ca[i], c[i], n[i+1]])
    # i = 1
    # rama.addBond([c[i-1], n[i], ca[i], c[i], n[i+1]])
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
    system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi, res_names = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi, oa.res_names

    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    parameter_i = []
    for i in range(nres):
        parameter_i.append(se_map_3_letter[res_names[i]])

    def lambda_coefficient(i, j, lambda_index):
        lambda_2_extra_terms = -0.5*alpha_coefficient(parameter_i[i],parameter_i[j],1)*np.log(p_antihb[parameter_i[i], parameter_i[j]][0])-0.25*alpha_coefficient(parameter_i[i], parameter_i[j], 2)*(np.log(p_antinhb[parameter_i[i+1],parameter_i[j-1]][0])+np.log(p_antinhb[parameter_i[i-1],parameter_i[j+1]][0]))-alpha_coefficient(parameter_i[i], parameter_i[j], 3)*(np.log(p_anti[parameter_i[i]])+np.log(p_anti[parameter_i[j]]))
        lambda_3_extra_terms = -alpha_coefficient(parameter_i[i],parameter_i[j], 4)*np.log(p_parhb[parameter_i[i+1],parameter_i[j]][0])-alpha_coefficient(parameter_i[i],parameter_i[j],5)*np.log(p_par[parameter_i[i+1]])+alpha_coefficient(parameter_i[i],parameter_i[j],4)*np.log(p_par[parameter_i[j]])
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
        return 0.0
    def alpha_coefficient(i,j, alpha_index):
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
        return 0.0

    # add beta potential
    # setup parameters
    k_beta = 1.0
    lambda_1 = [0]*nres*nres
    lambda_2 = [0]*nres*nres
    lambda_3 = [0]*nres*nres
    for i in range(1,nres-1):
        for j in range(1,nres-1):
            lambda_1[i+j*nres] = 1 #lambda_coefficient(i,j,1)
            lambda_2[i+j*nres] = 1 #lambda_coefficient(i,j,2)
            lambda_3[i+j*nres] = 1 #lambda_coefficient(i,j,3)

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

def read_beta_parameters(parameter_directory='.'):
    os.chdir(parameter_directory)
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

def apply_pap_term(oa):
    system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi, res_names = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi, oa.res_names

    pap_function = "-k_pap*gamma*0.5*(1+tanh(eta_pap*(r0-distance(p1,p2))))*0.5*(1+tanh(eta_pap*(r0-distance(p3,p4))))"
    # setup parameters
    k_pap = 10.0
    r0 = 0.8 # nm
    eta_pap = 70 # nm^-1
    gamma_aph = 1.0
    gamma_ap = 0.4
    gamma_p = 0.4

    pap = CustomCompoundBondForce(4, pap_function)
    pap.addGlobalParameter("k_pap", k_pap)
    pap.addGlobalParameter("r0", r0)
    pap.addGlobalParameter("eta_pap", eta_pap)
    pap.addPerBondParameter("gamma")

    for i in range(nres):
        for j in range(nres):
            # anti-parallel hairpin for i from 1 to N-13 and j from i+13 to min(i+16,N)
            # CAi CAj CAi+4 CAj-4
            # 1   2   3     4
            if i <= nres-13 and j >= i+13 and j <= min(i+16,nres):
                pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_aph])
            # anti-parallel for i from 1 to N-17 and j from i+17 to N
            # CAi CAj CAi+4 CAj-4
            # 1   2   3     4
            if i <= nres-17 and j >= i+17 and j <= nres:
                pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_ap])
            # parallel for i from 1 to N-13 and j from i+9 to N-4
            # CAi CAj CAi+4 CAj+4
            # 1   2   3     4
            if i <= nres-13 and j >= i+9 and j <= nres-4:
                pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_p])

    system.addForce(pap)

def apply_dsb_term(oa):
    system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi, res_names = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi, oa.res_names

    k_dsb = 1.0
    dsb_cutoff = .15
    eta_dsb = 100
    r_min_dsb = 0.6
    r_max_dsb = 0.7
    shift = {'ALA': 0.00*10, 'ARG': 2.04*10, 'ASN': 0.57*10, 'ASP': 0.57*10, 'CYS': 0.36*10, 'GLN': 1.11*10, 'GLU': 1.17*10, 'GLY': -1.52*10,   'HIS': 0.87*10, 'ILE': 0.67*10, 'LEU': 0.79*10, 'LYS': 1.47*10, 'MET': 1.03*10, 'PHE': 1.00*10, 'PRO': 0.10*10, 'SER': 0.26*10,  'THR': 0.37*10, 'TRP': 1.21*10, 'TYR': 1.15*10, 'VAL': 0.39*10}

    dsb = CustomNonbondedForce("k_dsb*0.5*(tanh(eta_dsb*(r-(r_min+shift1+shift2)))+tanh(eta_dsb*((r_max+shift1+shift2)-r)))")
    dsb.addGlobalParameter("k_dsb", k_dsb)
    dsb.addGlobalParameter("eta_dsb", eta_dsb)
    dsb.addGlobalParameter("r_min", r_min_dsb)
    dsb.addGlobalParameter("r_max", r_max_dsb)
    dsb.addPerParticleParameter("shift")
    for i in range(natoms):
        dsb.addParticle([shift[res_names[resi[i]]]])
    cb_with_gly_ca = [x if x >= 0 else y for x,y in zip(cb,ca)]
    dsb.addInteractionGroup(cb_with_gly_ca, cb_with_gly_ca)
    dsb.setCutoffDistance(dsb_cutoff)

    # find all pairs to include in the force
    interactions = []
    for cbi, cbj in combinations(cb_with_gly_ca, 2):
        if abs(resi[cbi]-resi[cbj]) >= 10:
            interactions.append((cbi, cbj))

    # get list of exclusions
    exclusions = get_exclusions(oa, interactions)

    # apply exclusions
    dsb.createExclusionsFromBonds(exclusions, 1)

    system.addForce(dsb)

def apply_helix_term(oa):
    system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi, res_names = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi, oa.res_names
    helix_frequencies = {
    'ALA': 0.77,
    'ARG': 0.68,
    'ASN': 0.07,
    'ASP': 0.15,
    'CYS': 0.23,
    'GLN': 0.33,
    'GLU': 0.27,
    'GLY': 0.0,
    'HIS': 0.06,
    'ILE': 0.23,
    'LEU': 0.62,
    'LYS': 0.65,
    'MET': 0.50,
    'PHE': 0.41,
    'PRO': 0.4,
    'SER': 0.35,
    'THR': 0.11,
    'TRP': 0.45,
    'TYR': 0.17,
    'VAL': 0.14
    }

    k_helix = 0
    gamma_prot = 2.0
    gamma_wat = -1.0
    eta_sigma = 70
    r_ON = .298
    r_OH = .206
    sigma_ON = .068
    sigma_OH = .076
    eta_direct = 50
    r_density_min = 0.45
    r_density_max = 0.65
    rho0 = 2.6

    rhoi = ""
    rhoip4 = ""
    for j in range(1,nres+1-3):
        rhoi += "0.25*(1+tanh(eta_direct*(distance(p4,p%d)-r_density_min)))*(1+tanh(eta_direct*(r_density_max-distance(p4,p%d))))+" % (j+5, j+5)
        rhoip4 += "0.25*(1+tanh(eta_direct*(distance(p5,p%d)-r_density_min)))*(1+tanh(eta_direct*(r_density_max-distance(p5,p%d))))+" % (j+5, j+5)
    rhoi = rhoi[:-1]
    rhoip4 = rhoip4[:-1]

    helix_function = "-k_helix*(fai+faip4)*(gamma_prot*(1-sigma_wat)+gamma_wat*sigma_wat)*exp(-(distance(p1,p2)-r_ON)^2/2*sigma_ON^2-(distance(p1,p3)-r_OH)^2/2*sigma_OH^2);sigma_wat=0.25*(1-tanh(eta_sigma*(rhoi-rho0))*(1-tanh(eta_sigma*(rhoip4-rho0))));rhoi=%s;rhoip4=%s" % (rhoi, rhoip4)

    helix = CustomCompoundBondForce(5+nres-3, helix_function)
    helix.addGlobalParameter("k_helix", k_helix)
    helix.addGlobalParameter("gamma_prot", gamma_prot)
    helix.addGlobalParameter("gamma_wat", gamma_wat)
    helix.addGlobalParameter("eta_sigma", eta_sigma)
    helix.addGlobalParameter("eta_direct", eta_direct)
    helix.addGlobalParameter("sigma_ON", sigma_ON)
    helix.addGlobalParameter("sigma_OH", sigma_OH)
    helix.addGlobalParameter("r_density_min", r_density_min)
    helix.addGlobalParameter("r_density_max", r_density_max)
    helix.addGlobalParameter("rho0", rho0)
    helix.addGlobalParameter("r_ON", r_ON)
    helix.addGlobalParameter("r_OH", r_OH)
    helix.addPerBondParameter("i")
    helix.addPerBondParameter("ip4")
    helix.addPerBondParameter("fai")
    helix.addPerBondParameter("faip4")

    cb_with_gly_ca = [x if x >= 0 else y for x,y in zip(cb,ca)]
    #for i in range(nres-4):
    for i in range(2,3):
        print(i)
        other_cb_with_gly_ca = list(cb_with_gly_ca)
        print(other_cb_with_gly_ca)
        del other_cb_with_gly_ca[i-1]
        print(other_cb_with_gly_ca)
        del other_cb_with_gly_ca[i-1]
        print(other_cb_with_gly_ca)
        del other_cb_with_gly_ca[i-1]
        print(other_cb_with_gly_ca)
        if res_names[i+4] == "PRO":
            continue
        fai = helix_frequencies[res_names[i]]
        faip4 = helix_frequencies[res_names[i+4]]
        helix.addBond([o[i], n[i+4], h[i+4], cb_with_gly_ca[i], cb_with_gly_ca[i+4], *other_cb_with_gly_ca], [i+1, i+5, fai, faip4])

    system.addForce(helix)

def read_memory(pdb_file, chain_name, target_start, fragment_start, length, weight, min_seq_sep, max_seq_sep, ca, cb):
    memory_interactions = []

    if not os.path.isfile(pdb_file):
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(pdb_file.split('.')[0].lower(), pdir='.')
        os.rename("pdb%s.ent" % pdb_id, "%s.pdb" % pdb_id)

    p = PDBParser()
    structure = p.get_structure('X', pdb_file)
    chain = structure[0][chain_name]
    residues = [x for x in chain if x.get_full_id()[3][1] in range(fragment_start,fragment_start+length-1)]
    for i, residue_i in enumerate(residues):
        for j, residue_j in enumerate(residues):
            if abs(i-j) > max_seq_sep:
                continue
            target_index_i = target_start + i - 1
            target_index_j = target_start + j - 1
            atom_list = []
            target_atom_list = []
            if abs(i-j) >= min_seq_sep:
                ca_i = residue_i['CA']
                atom_list.append(ca_i)
                target_atom_list.append(ca[target_index_i])
                ca_j = residue_j['CA']
                atom_list.append(ca_j)
                target_atom_list.append(ca[target_index_j])
                if not residue_i.get_resname() == "GLY" and cb[target_index_i] >= 0:
                    cb_i = residue_i['CB']
                    atom_list.append(cb_i)
                    target_atom_list.append(cb[target_index_i])
                if not residue_j.get_resname() == "GLY" and cb[target_index_j] >= 0:
                    cb_j = residue_j['CB']
                    atom_list.append(cb_j)
                    target_atom_list.append(cb[target_index_j])
            for atom_i, atom_j in combinations(atom_list, 2):
                particle_1 = target_atom_list[atom_list.index(atom_i)]
                particle_2 = target_atom_list[atom_list.index(atom_j)]
                r_ijm = (atom_i - atom_j)/10.0 # convert to nm
                sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = 1.0
                w_m = weight
                memory_interaction = [particle_1, particle_2, [w_m, gamma_ij, r_ijm, sigma_ij]]
                memory_interactions.append(memory_interaction)
    return memory_interactions

def apply_associative_memory_term(oa):
    system, nres, n, h, ca, c, o, cb, res_type, natoms, bonds, resi, res_names = oa.system, oa.nres, oa.n, oa.h, oa.ca, oa.c, oa.o, oa.cb, oa.res_type, oa.natoms, oa.bonds, oa.resi, oa.res_names

    k_am = 1.0
    min_seq_sep = 3
    max_seq_sep = 9
         #pdbid #chain #target #fragment #length #weight
    memories = [['9-peptide.pdb', 'C', 1, 1, 63, 1]]

    am_function = '-k_am*w_m*gamma_ij*exp(-(r-r_ijm)^2/(2*sigma_ij^2))'
    am = CustomBondForce(am_function)

    am.addGlobalParameter('k_am', k_am)
    am.addPerBondParameter('w_m')
    am.addPerBondParameter('gamma_ij')
    am.addPerBondParameter('r_ijm')
    am.addPerBondParameter('sigma_ij')

    for memory in memories:
        memory_interactions = read_memory(*memory, min_seq_sep, max_seq_sep, ca, cb)
        for memory_interaction in memory_interactions:
            am.addBond(*memory_interaction)

    system.addForce(am)

