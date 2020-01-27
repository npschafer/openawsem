from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from Bio.PDB.PDBParser import PDBParser

def read_reference_structure_for_q_calculation_3(oa, pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, Qflag=0, a=0.1):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
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

def qbias_term(oa, q0, reference_pdb_file, reference_chain_name, k_qbias=100*kilocalorie_per_mole, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.8*nanometers):
    k_qbias = k_qbias.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    qbias = CustomCVForce("0.5*k_qbias*(q-q0)^2")
    q = q_value(oa, reference_pdb_file, reference_chain_name, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
    qbias.addCollectiveVariable("q", q)
    qbias.addGlobalParameter("k_qbias", k_qbias)
    qbias.addGlobalParameter("q0", q0)
    return qbias
