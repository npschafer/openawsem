from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
import pandas as pd

def read_reference_structure_for_q_calculation_4(oa, contact_threshold,rnative_dat,  min_seq_sep=3, max_seq_sep=np.inf):
    # use contact matrix for Q calculation
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
    in_rnative = np.loadtxt(rnative_dat)  ## read in rnative_dat file for Q calculation
    structure_interactions = []
    chain_start = 0
    count = 0;
    for i in range(oa.nres):
        chain_start += count
        count = 0
        for j in range(oa.nres):
            count +=1
            if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                r_ijN = in_rnative[i][j]/10.0 * nanometers # convert to nm
                if r_ijN < contact_threshold: continue
                sigma_ij = 0.1*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = 1.0
                i_index = oa.ca[i]
                j_index = oa.ca[j]
                structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                structure_interactions.append(structure_interaction)
    return structure_interactions



def q_value_dat(oa, contact_threshold, rnative_dat="rnative.dat", min_seq_sep=3, max_seq_sep=np.inf):
    ### Added by Mingchen
    ### this function is solely used for template based modelling from rnative.dat file
    ### for details, refer to Chen, Lin & Lu Wolynes JCTC 2018
    qvalue_dat = CustomBondForce("(1/normalization)*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue_dat.addPerBondParameter("gamma_ij")
    qvalue_dat.addPerBondParameter("r_ijN")
    qvalue_dat.addPerBondParameter("sigma_ij")
    structure_interactions_tbm_q = read_reference_structure_for_q_calculation_4(oa, contact_threshold=contact_threshold,rnative_dat="rnative.dat",  min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep)
    qvalue_dat.addGlobalParameter("normalization", len(structure_interactions_tbm_q))
    for structure_interaction_tbm_q in structure_interactions_tbm_q:
        qvalue_dat.addBond(*structure_interaction_tbm_q)
    return qvalue_dat


def tbm_q_term(oa, k_tbm_q, tbm_q_min_seq_sep=2, tbm_q_cutoff=0.2*nanometers, tbm_q_well_width=0.1, target_q=1.0):
    ### Added by Mingchen Chen
    ### this function is solely used for template based modelling from rnative.dat file
    ### for details, refer to Chen, Lin & Lu Wolynes JCTC 2018
    print("TBM_Q term ON")
    tbm_q = CustomCVForce("0.5*k_tbm_q*(q-q0)^2")
    q = q_value_dat(oa, contact_threshold=tbm_q_cutoff, rnative_dat="rnative.dat", min_seq_sep=tbm_q_min_seq_sep, max_seq_sep=np.inf)
    tbm_q.addCollectiveVariable("q", q)
    tbm_q.addGlobalParameter("k_tbm_q", k_tbm_q)
    tbm_q.addGlobalParameter("q0", target_q)
    tbm_q.setForceGroup(22)
    return tbm_q



def fragment_memory_term(oa, k_fm=0.04184, frag_file_list_file="./frag.mem", npy_frag_table="./frag_table.npy",
                    min_seq_sep=3, max_seq_sep=9, fm_well_width=0.1, UseSavedFragTable=True):
    # 0.8368 = 0.01 * 4.184 # in kJ/mol, converted from default value in LAMMPS AWSEM
    k_fm *= oa.k_awsem
    frag_table_rmin = 0
    frag_table_rmax = 5  # in nm
    frag_table_dr = 0.01
    r_array = np.arange(frag_table_rmin, frag_table_rmax, frag_table_dr)
    number_of_atoms = oa.natoms
    r_table_size = int((frag_table_rmax - frag_table_rmin)/frag_table_dr)  # 500 here.
    raw_frag_table = np.zeros((number_of_atoms, 6*(1+max_seq_sep), r_table_size))
    data_dic = {}
    for i in range(oa.natoms):
        if i in oa.ca:
            res_id = oa.resi[i]    # oa.resi start with 0, but pdb residue id start with 1
            data_dic[("CA", 1+int(res_id))] = i
        if i in oa.cb:
            res_id = oa.resi[i]
            data_dic[("CB", 1+int(res_id))] = i
    # print(oa.res_type)
    # print(oa.resi)
    # print(data_dic)
    frag_location_pre = os.path.dirname(frag_file_list_file)
    # frag_file_list_file = frag_location_pre + "frags.mem"
    # frag_table_file = frag_location_pre + "frag_table.npy"
    frag_table_file = npy_frag_table

    if os.path.isfile(frag_table_file) and UseSavedFragTable:
        print(f"Reading Fragment table from {frag_table_file}.")
        frag_table, interaction_list, interaction_pair_to_bond_index = np.load(frag_table_file, allow_pickle=True)
        print(f"Fragment table loaded, number of bonds: {len(interaction_list)}")
        frag_file_list = []
    else:
        print(f"Fragment table file is not found. Reading fragments files.")
        frag_file_list = pd.read_csv(frag_file_list_file, skiprows=4, sep="\s+", names=["location", "target_start", "fragment_start", "frag_len", "weight"])
        interaction_list = set()
    for frag_index in range(len(frag_file_list)):
        location = frag_file_list["location"].iloc[frag_index]
        frag_name = os.path.join(frag_location_pre, location)
        frag_len = frag_file_list["frag_len"].iloc[frag_index]
        weight = frag_file_list["weight"].iloc[frag_index]
        target_start = frag_file_list["target_start"].iloc[frag_index]  # residue id
        fragment_start = frag_file_list["fragment_start"].iloc[frag_index]  # residue id
        frag = pd.read_csv(frag_name, skiprows=2, sep="\s+", header=None, names=["Res_id", "Res", "Type", "i", "x", "y", "z"])
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
                if seq_sep > max_seq_sep:
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
    if (not os.path.isfile(frag_table_file)) or (not UseSavedFragTable):
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
    # for i in range(oa.natoms):
    #     fm.addParticle([i])

    # # add interaction that are cutoff away
    # # print(sorted(interaction_list))
    # for (i, j) in interaction_list:
    #     fm.addInteractionGroup([i], [j])
    # # add per-particle parameters
    # fm.addPerParticleParameter("index")

    # for edge case, that r > frag_table_rmax
    max_r_index_1 = r_table_size - 2
    fm = CustomCompoundBondForce(2, f"-k_fm*((v2-v1)*r+v1*r_2-v2*r_1)/(r_2-r_1); \
                                v1=frag_table(index, r_index_1);\
                                v2=frag_table(index, r_index_2);\
                                r_1=frag_table_rmin+frag_table_dr*r_index_1;\
                                r_2=frag_table_rmin+frag_table_dr*r_index_2;\
                                r_index_2=r_index_1+1;\
                                r_index_1=min({max_r_index_1}, floor(r/frag_table_dr));\
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


def read_memory(oa, pdb_file, chain_name, target_start, fragment_start, length, weight, min_seq_sep, max_seq_sep, am_well_width=0.1):
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
                target_atom_list_i.append(oa.ca[target_index_i])
                ca_j = residue_j['CA']
                atom_list_j.append(ca_j)
                target_atom_list_j.append(oa.ca[target_index_j])
                if not residue_i.get_resname() == "GLY" and oa.cb[target_index_i] >= 0:
                    cb_i = residue_i['CB']
                    atom_list_i.append(cb_i)
                    target_atom_list_i.append(oa.cb[target_index_i])
                if not residue_j.get_resname() == "GLY" and oa.cb[target_index_j] >= 0:
                    cb_j = residue_j['CB']
                    atom_list_j.append(cb_j)
                    target_atom_list_j.append(oa.cb[target_index_j])
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

def associative_memory_term(oa, memories, k_am=0.8368, min_seq_sep=3, max_seq_sep=9, am_well_width=0.1):
    # 0.8368 = 0.2 * 4.184 # in kJ/mol, converted from default value in LAMMPS AWSEM
    #pdbid #chain #target #fragment #length #weight
    # multiply interaction strength by overall scaling
    k_am *= oa.k_awsem
    am_function = '-k_am*w_m*gamma_ij*exp(-(r-r_ijm)^2/(2*sigma_ij^2))'
    am = CustomBondForce(am_function)
    am.addGlobalParameter('k_am', k_am)
    am.addPerBondParameter('w_m')
    am.addPerBondParameter('gamma_ij')
    am.addPerBondParameter('r_ijm')
    am.addPerBondParameter('sigma_ij')
    for memory in memories:
        memory_interactions = read_memory(oa, *memory, min_seq_sep, max_seq_sep, am_well_width=am_well_width)
        for memory_interaction in memory_interactions:
            am.addBond(*memory_interaction)
    return am



def density_dependent_associative_memory_term(oa, memories, k_am_dd=1.0, am_dd_min_seq_sep=3, am_dd_max_seq_sep=9, eta_density=50, r_density_min=.45, r_density_max=.65, density_alpha=1.0, density_normalization=2.0, rho0=2.6, am_well_width=0.1, density_min_seq_sep=10, density_only_from_native_contacts=False, density_pdb_file=None, density_chain_name=None, density_native_contact_min_seq_sep=4, density_native_contact_threshold=0.8*nanometers):

    k_am_dd *= oa.k_awsem

    am_dd = CustomGBForce()

    # add all particles to force
    for i in range(oa.natoms):
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
        structure_interactions = read_amhgo_structure(oa, pdb_file=density_pdb_file, chain_name=density_chain_name, amhgo_min_seq_sep=density_native_contact_min_seq_sep, amhgo_contact_threshold=density_native_contact_threshold, amhgo_well_width=0.1) # the well width is not used, so the value doesn't matter

        native_contacts = []
        for interaction in structure_interactions:
            i_index, j_index, [gamma_ij, r_ijN, sigma_ij] = interaction
            native_contacts.append((i_index, j_index))
            native_contacts.append((j_index, i_index))

    # setup tabulated functions and interactions
    density_gamma_ij = [0.0]*oa.natoms*oa.natoms
    for i in range(oa.natoms):
        for j in range(oa.natoms):
            if (i in oa.cb or (oa.res_type[oa.resi[i]] == "IGL" and i in oa.ca)) and (j in oa.cb or (oa.res_type[oa.resi[j]] == "IGL" and i in oa.ca)) and abs(oa.resi[i]-oa.resi[j])>=density_min_seq_sep:
                if not density_only_from_native_contacts or (i, j) in native_contacts or (j, i) in native_contacts:
                    density_gamma_ij[i+j*oa.natoms] = 1.0
                    density_gamma_ij[j+i*oa.natoms] = 1.0
    am_dd.addTabulatedFunction("density_gamma_ij", Discrete2DFunction(oa.natoms, oa.natoms, density_gamma_ij))

    gamma_ij = [0.0]*oa.natoms*oa.natoms*len(memories)
    sigma_ij = [0.1]*oa.natoms*oa.natoms*len(memories)
    r_ijm = [0.0]*oa.natoms*oa.natoms*len(memories)
    for k, memory in enumerate(memories):
        memory_interactions = read_memory(oa, *memory, am_dd_min_seq_sep, am_dd_max_seq_sep, am_well_width=am_well_width)
        for memory_interaction in memory_interactions:
            i, j, (w_m, gamma, r, sigma) = memory_interaction
            gamma_ij[i+j*oa.natoms+k*oa.natoms*oa.natoms] = gamma
            gamma_ij[j+i*oa.natoms+k*oa.natoms*oa.natoms] = gamma
            sigma_ij[i+j*oa.natoms+k*oa.natoms*oa.natoms] = sigma
            sigma_ij[j+i*oa.natoms+k*oa.natoms*oa.natoms] = sigma
            r_ijm[i+j*oa.natoms+k*oa.natoms*oa.natoms] = r
            r_ijm[j+i*oa.natoms+k*oa.natoms*oa.natoms] = r
    am_dd.addTabulatedFunction("gamma_ij", Discrete3DFunction(oa.natoms, oa.natoms, len(memories), gamma_ij))
    am_dd.addTabulatedFunction("sigma_ij", Discrete3DFunction(oa.natoms, oa.natoms, len(memories), sigma_ij))
    am_dd.addTabulatedFunction("r_ijm", Discrete3DFunction(oa.natoms, oa.natoms, len(memories), r_ijm))

    # add computed values
    # compute the density
    am_dd.addComputedValue("rho", "0.25*density_gamma_ij(index1, index2)*(1+tanh(eta_density*(r-r_density_min)))*(1+tanh(eta_density*(r_density_max-r)))", CustomGBForce.ParticlePair)

    # function that determines how the AM term depends on density
    #f_string = "0.25*(1-tanh(eta_density*(rho0-rho1)))*(1-tanh(eta_density*(rho0-rho2)))" # both residues must be buried for the interaction to be active
    f_string = "1-(0.25*(1-tanh(eta_density*(rho1-rho0)))*(1-tanh(eta_density*(rho2-rho0))))" # one residue being buried is enough for the interaction to be active

    # add energy term for each memory
    for k, memory in enumerate(memories):
        memory_interactions = read_memory(oa, *memory, am_dd_min_seq_sep, am_dd_max_seq_sep, am_well_width=am_well_width)
        for memory_interaction in memory_interactions:
            i, j, (w_m, gamma, r, sigma) = memory_interaction
        am_dd.addEnergyTerm("-k_am_dd*(density_alpha*f*density_normalization*beta_ij+(1-density_alpha)*beta_ij);\
        beta_ij=%f*gamma_ij(index1,index2,%d)*exp(-(r-r_ijm(index1,index2,%d))^2/(2*sigma_ij(index1,index2,%d)^2));\
        f=%s" % (w_m, k, k, k, f_string), CustomGBForce.ParticlePair)

    return am_dd

def read_amhgo_structure(oa, pdb_file, chain_name, amhgo_min_seq_sep=4, amhgo_contact_threshold=0.8*nanometers, amhgo_well_width=0.1):
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
                            i_index = oa.ca[i]
                        if atom_i in cb_list:
                            i_index = oa.cb[i]
                        if atom_j in ca_list:
                            j_index = oa.ca[j]
                        if atom_j in cb_list:
                            j_index = oa.cb[j]
                        structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                        print(i_index, j_index, gamma_ij, r_ijN, sigma_ij)
                        structure_interactions.append(structure_interaction)
    return structure_interactions

def additive_amhgo_term(oa, pdb_file, chain_name, k_amhgo=4.184, amhgo_min_seq_sep=10, amhgo_contact_threshold=0.8*nanometers, amhgo_well_width=0.1):
    import itertools
    # multiply interaction strength by overall scaling
    print("AMH-GO structure based term is ON")
    k_amhgo *= oa.k_awsem
    # create contact force
    amhgo = CustomBondForce("-k_amhgo*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    # # add global parameters
    amhgo.addGlobalParameter("k_amhgo", k_amhgo)
    amhgo.addPerBondParameter("gamma_ij")
    amhgo.addPerBondParameter("r_ijN")
    amhgo.addPerBondParameter("sigma_ij")
    # create bonds
    structure_interactions = read_amhgo_structure(oa, pdb_file, chain_name, amhgo_min_seq_sep, amhgo_contact_threshold, amhgo_well_width=amhgo_well_width)
    print(structure_interactions)
    for structure_interaction in structure_interactions:
        print(structure_interaction)
        amhgo.addBond(*structure_interaction)
    #amhgo.setForceGroup(22)
    return amhgo

def er_term(oa, k_er=4.184, er_min_seq_sep=2, er_cutoff=99.0, er_well_width=0.1):
    ### this is a structure prediction related term; Adapted from Sirovitz Schafer Wolynes 2017 Protein Science;
    ### See original papers for reference: Make AWSEM AWSEM-ER with Evolutionary restrictions
    ### ER restrictions can be obtained from multiple sources (RaptorX, deepcontact, and Gremlin)
    ### term modified from amh-go term, and the current strength seems to be high, and needs to be lowered somehow.
    ### amh-go normalization factor will be added soon. Based on Eastwood Wolynes 2000 JCP
    print("ER term is ON")
    import itertools
    k_er *= oa.k_awsem
    # create contact force
    er = CustomBondForce("-k_er*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    # # add global parameters
    er.addGlobalParameter("k_er", k_er)
    er.addPerBondParameter("gamma_ij")
    er.addPerBondParameter("r_ijN")
    er.addPerBondParameter("sigma_ij")
    structure_interactions_er = []
    ### read in dat files from contact predictions;
    in_rnativeCACA = np.loadtxt('go_rnativeCACA.dat')
    in_rnativeCACB = np.loadtxt('go_rnativeCACB.dat')
    in_rnativeCBCB = np.loadtxt('go_rnativeCBCB.dat')
    for i in range(oa.nres):
        for j in range(oa.nres):
            if abs(i-j) >= er_min_seq_sep and in_rnativeCACA[i][j]<er_cutoff:
                sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = 1.0
                r_ijN = in_rnativeCACA[i][j]/10.0*nanometers;
                structure_interactions_er.append([oa.ca[i], oa.ca[j], [gamma_ij, r_ijN, sigma_ij]])
            if abs(i-j) >= er_min_seq_sep and in_rnativeCACB[i][j]<er_cutoff and oa.cb[j]!= -1:
                sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = 1.0
                r_ijN = in_rnativeCACB[i][j]/10.0*nanometers;
                structure_interactions_er.append([oa.ca[i], oa.cb[j], [gamma_ij, r_ijN, sigma_ij]])
            if abs(i-j) >= er_min_seq_sep and in_rnativeCBCB[i][j]<er_cutoff and oa.cb[j]!= -1 and oa.cb[i]!= -1:#oa.res_type[oa.resi[i]] != "IGL" and oa.res_type[oa.resi[j]] != "IGL":
                sigma_ij = er_well_width*abs(i-j)**0.15 # 0.1 nm = 1 A
                gamma_ij = 1.0
                r_ijN = in_rnativeCBCB[i][j]/10.0*nanometers;
                structure_interactions_er.append([oa.cb[i], oa.cb[j], [gamma_ij, r_ijN, sigma_ij]])
                # print([i, j, oa.res_type[oa.resi[i]], oa.res_type[oa.resi[j]],oa.cb[i], oa.cb[j], [gamma_ij, r_ijN, sigma_ij]])
    # create bonds
    for structure_interaction_er in structure_interactions_er:
        er.addBond(*structure_interaction_er)
    er.setForceGroup(21)
    return er



'''
# will be deleted in the future.
def read_reference_structure_for_q_calculation(oa, pdb_file, chain_name, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.8*nanometers):
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
                            i_index = oa.ca[i]
                        if atom_i in cb_list:
                            i_index = oa.cb[i]
                        if atom_j in ca_list:
                            j_index = oa.ca[j]
                        if atom_j in cb_list:
                            j_index = oa.cb[j]
                        structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                        structure_interactions.append(structure_interaction)

    return structure_interactions

def read_reference_structure_for_q_calculation_2(oa, pdb_file, min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.8*nanometers):
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
                                i_index = oa.ca[i+chain_start]
                            if atom_i in cb_list:
                                i_index = oa.cb[i+chain_start]
                            if atom_j in ca_list:
                                j_index = oa.ca[j+chain_start]
                            if atom_j in cb_list:
                                j_index = oa.cb[j+chain_start]
                            structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                            structure_interactions.append(structure_interaction)

    return structure_interactions
'''