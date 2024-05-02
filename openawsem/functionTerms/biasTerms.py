try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
import numpy as np
from Bio.PDB.PDBParser import PDBParser

def read_reference_structure_for_q_calculation_3(oa, pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, Qflag=0, a=0.1, removeDNAchains=True):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
    structure_interactions = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    model = structure[0]
    chain_start = 0
    count = 0
    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    proteinResidues += ["NGP", "IGL", "IPR"]
    rnaResidues = ['A', 'G', 'C', 'U', 'I']
    dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']

    for chain in model.get_chains():
        chain_start += count
        count = 0
        if removeDNAchains and np.alltrue([a.get_resname().strip() in dnaResidues for a in chain.get_residues()]):
            print(f"chain {chain.id} is a DNA chain. will be ignored for Q evaluation")
            continue
        elif removeDNAchains and np.alltrue([a.get_resname().strip() not in proteinResidues for a in chain.get_residues()]):
            print(f"chain {chain.id} is a ligand chain. will be ignored for Q evaluation")
            continue
        # print(chain)
        for i, residue_i in enumerate(chain.get_residues()):
            #  print(i, residue_i)
            count +=1
            for j, residue_j in enumerate(chain.get_residues()):
                if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                    ca_i = residue_i['CA']

                    ca_j = residue_j['CA']

                    r_ijN = abs(ca_i - ca_j)/10.0*nanometers # convert to nm
                    if Qflag ==1 and r_ijN >= contact_threshold: continue
                    sigma_ij = a*(abs(i-j)**0.15)  # 0.1 nm = 1 A
                    gamma_ij = 1.0

                    if reference_chain_name != "ALL" and (chain.id not in reference_chain_name):
                        continue
                    i_index = oa.ca[i+chain_start]
                    j_index = oa.ca[j+chain_start]
                    structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                    structure_interactions.append(structure_interaction)
    # print("Done reading")
    # print(structure_interactions)
    return structure_interactions

def read_reference_structure_for_qc_calculation(oa, pdb_file, min_seq_sep=3, a=0.1, startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=None):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
    structure_interactions = []
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    if endResidueIndex == -1:
        endResidueIndex = len(list(structure.get_residues()))
    if residueIndexGroup is None:
        # residueIndexGroup is used for non-continuous residues that used for Q computation.
        residueIndexGroup = list(range(len(list(structure.get_residues()))))
    for i, res_i in enumerate(structure.get_residues()):
        chain_i = res_i.get_parent().id
        if i < startResidueIndex:
            continue
        if i not in residueIndexGroup:
            continue
        for j, res_j in enumerate(structure.get_residues()):
            if j > endResidueIndex:
                continue
            if j not in residueIndexGroup:
                continue
            chain_j = res_j.get_parent().id
            if j-i >= min_seq_sep and chain_i == chain_j:
                    ca_i = res_i['CA']
                    ca_j = res_j['CA']

                    r_ijN = abs(ca_i - ca_j)/10.0  # convert to nm
                    sigma_ij = a*(abs(i-j)**0.15)  # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    i_index = oa.ca[i]
                    j_index = oa.ca[j]
                    structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                    structure_interactions.append(structure_interaction)
    return structure_interactions

def q_value(oa, reference_pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, forceGroup=1):
    ### Modified by Mingchen to compute canonical QW/QO

    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_q_calculation_3(oa, reference_pdb_file, reference_chain_name=reference_chain_name,
        min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=0)
    # print(len(structure_interactions))
    # print(structure_interactions)
    # create bond force for q calculation
    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")

    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(forceGroup)
    return qvalue


def qc_value(oa, reference_pdb_file, min_seq_sep=10, a=0.2):
    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_qc_calculation(oa, reference_pdb_file, min_seq_sep=min_seq_sep, a=a)
    # print(len(structure_interactions))
    # print(structure_interactions)

    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(3)
    return qvalue

def partial_q_value(oa, reference_pdb_file, min_seq_sep=3, a=0.1, startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=None, forceGroup=4):
    print(f"Including partial q value computation, start residue index: {startResidueIndex}, end residue index: {endResidueIndex}, residueIndexGroup: {residueIndexGroup}")
    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_qc_calculation(oa, reference_pdb_file, min_seq_sep=min_seq_sep, a=a, startResidueIndex=startResidueIndex, endResidueIndex=endResidueIndex, residueIndexGroup=residueIndexGroup)
    # print(len(structure_interactions))
    # print(structure_interactions)
    if len(structure_interactions) == 0:
        print("No atom found, Please check your startResidueIndex and endResidueIndex.")
        exit()
    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(forceGroup)
    return qvalue

def qbias_term(oa,reference_pdb_file, q0, reference_chain_name="ALL", k_qbias=200*kilocalorie_per_mole, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.8*nanometers, forceGroup=4):
    k_qbias = k_qbias.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    qbias = CustomCVForce(f"0.5*{k_qbias}*(q-{q0})^2")
    # qbias = CustomCVForce(f"0.5*{k_qbias}*(q-q0)^2")
    q = q_value(oa, reference_pdb_file, reference_chain_name, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
    if oa.periodic:
        q.setUsesPeriodicBoundaryConditions(True)
    qbias.addCollectiveVariable("q", q)
    # qbias.addGlobalParameter("k_qbias", k_qbias)
    # qbias.addGlobalParameter("q0", q0)
    is_periodic=qbias.usesPeriodicBoundaryConditions()
    print("\nqbias_term is in PBC",is_periodic)
    print("\nqbias is",q0)
    qbias.setForceGroup(forceGroup)
    return qbias

def create_dist(oa,fileA="groupA.dat",fileB="groupB.dat",forceGroup=4):
    #groupA = list(range(68))
    #groupB = list(range(196, oa.natoms))
    #print(cA)
    #print(cB)
    groupA = np.array(np.loadtxt(fileA,dtype=int)).tolist()
    #print (groupA)
    #groupA = [0,1]
    groupB = np.array(np.loadtxt(fileB,dtype=int)).tolist()
    #groupA, groupB = get_contact_atoms('crystal_structure-openmmawsem.pdb', chainA=cA, chainB=cB)
    #pull_d = CustomCentroidBondForce(2, 'distance(g1,g2)-R0') # 为什么这里给了R0，实际distance变成了2倍？
    #pull_d.addGlobalParameter("R0", 0.0*angstroms)
    pull_d = CustomCentroidBondForce(2, 'distance(g1,g2)')
    pull_d.addGroup(groupA)
    pull_d.addGroup(groupB) # addGroup(groupB)
    pull_d.addBond([0, 1])
    pull_d.setForceGroup(forceGroup)
    return pull_d

def create_centroid_system(oa, fileA="groupA.dat",fileB="groupB.dat",k=100,R0=0,forceGroup=26):
    K_pull = k*4.184 * oa.k_awsem
    #R0 = R0*nm
    pull_force = CustomCVForce("0.5*K_pull*(d-R0)^2")
    d = create_dist(oa)#    cA = list(range(68)),
    pull_force.addCollectiveVariable("d", d)
    pull_force.addGlobalParameter("K_pull", K_pull)
    pull_force.addGlobalParameter("R0", R0)
    pull_force.setForceGroup(forceGroup)
    return pull_force

def create_dist_vector(oa,fileA="groupA.dat",fileB="groupB.dat",fileC="groupC.dat",forceGroup=4):
    #groupA = list(range(68))
    #groupB = list(range(196, oa.natoms))
    #print(cA)
    #print(cB)
    groupA = np.array(np.loadtxt(fileA,dtype=int)).tolist()
    groupB = np.array(np.loadtxt(fileB,dtype=int)).tolist()
    groupC = np.array(np.loadtxt(fileC,dtype=int)).tolist()
    print (groupA)
    print (groupB)
    print (groupC)
    #groupA, groupB = get_contact_atoms('crystal_structure-openmmawsem.pdb', chainA=cA, chainB=cB)
    #pull_d = CustomCentroidBondForce(2, 'distance(g1,g2)-R0') # 为什么这里给了R0，实际distance变成了2倍？
    #pull_d.addGlobalParameter("R0", 0.0*angstroms)
    pull_d = CustomCentroidBondForce(3, f"r1*cos(theta);\
                                r1=distance(p1,p2);\
                                theta=angle(p1, p2, p3);")

    #pull_d = CustomCompoundBondForce(3, f"r1*cos(theta);\
     #                           r1=distance(p1,p2);\
      #                          theta=angle(p1, p2, p3);")

    g1=pull_d.addGroup(groupA)
    g2=pull_d.addGroup(groupB) # addGroup(groupB)
    #pull_d.addBond([0,1,2])
    g3=pull_d.addGroup(groupC)
    #g4=pull_d.addGroup(groupA)
    #g5=pull_d.addGroup(groupB)


    #print (g1,g2,g3,g4)
    #pull_d.addBond([2800,3435,3780])
    pull_d.addBond([0,1,2],[])
    #pull_d.addBond([g0])
    #pull_d.addBond([1, 2])
    #pull_d.addBond([0, 3])
    pull_d.setForceGroup(forceGroup)
    return pull_d

def create_centroid_system2(oa, fileA="groupA.dat",fileB="groupB.dat",fileC="groupC.dat",k=100,R0=0,forceGroup=26):
    K_pull = k*4.184 * oa.k_awsem
    #R0 = R0*nm
    pull_force = CustomCVForce("0.5*K_pull*(d-R0)^2")
    d = create_dist_vector(oa,fileA=fileA,fileB=fileB,fileC=fileC)#    cA = list(range(68)),
    pull_force.addCollectiveVariable("d", d)
    pull_force.addGlobalParameter("K_pull", K_pull)
    pull_force.addGlobalParameter("R0", R0)
    pull_force.setForceGroup(forceGroup)
    return pull_force

def rg_term(oa, convertToAngstrom=True):
    rg_square = CustomBondForce("1/normalization*r^2")
    # rg = CustomBondForce("1")
    rg_square.addGlobalParameter("normalization", oa.nres*oa.nres)
    for i in range(oa.nres):
        for j in range(i+1, oa.nres):
            rg_square.addBond(oa.ca[i], oa.ca[j], [])
    if convertToAngstrom:
        unit = 10
    else:
        unit = 1
    rg = CustomCVForce(f"{unit}*rg_square^0.5")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(2)
    return rg

def rg_bias_term(oa, k=1*kilocalorie_per_mole, rg0=0, atomGroup=-1, forceGroup=27):
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_rg = oa.k_awsem * k
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup     # atomGroup = [0, 1, 10, 12]  means include residue 1, 2, 11, 13.
    n = len(group)
    normalization = n*n
    rg_square = CustomBondForce(f"1.0/{normalization}*r^2")
    # rg = CustomBondForce("1")
    # rg_square.addGlobalParameter("normalization", n*n)
    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond(ca[i], ca[j], [])
    rg = CustomCVForce(f"{k_rg}*(rg_square^0.5-{rg0})^2")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(forceGroup)
    return rg

def cylindrical_rg_bias_term(oa, k=1*kilocalorie_per_mole, rg0=0, atomGroup=-1, forceGroup=27):
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_rg = oa.k_awsem * k
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup          # atomGroup = [0, 1, 10, 12]  means include residue 1, 2, 11, 13.
    n = len(group)
    normalization = n * n
    rg_square = CustomCompoundBondForce(2, f"1/{normalization}*((x1-x2)^2+(y1-y2)^2)")

    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond([ca[i], ca[j]], [])

    rg = CustomCVForce(f"{k_rg}*(rg_square^0.5-{rg0})^2")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(forceGroup)
    return rg

def pulling_term(oa, k_pulling=4.184, forceDirect="x", appliedToResidue=1, forceGroup=19):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k_pulling *= oa.k_awsem
    pulling = CustomExternalForce(f"(-{k_pulling})*({forceDirect})")
    for i in range(oa.natoms):
        if appliedToResidue == "LAST":
            appliedToResidue = oa.nres
        if appliedToResidue == "FIRST":
            appliedToResidue = 1
        if oa.resi[i] == (appliedToResidue-1):
            pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    pulling.setForceGroup(forceGroup)
    return pulling