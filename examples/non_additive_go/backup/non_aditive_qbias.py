from openmmawsem import *
def non_aditive_qbias_term(oa, 
                          q0,
                          reference_pdb_file, 
                          reference_chain_name, 
                          k_qbias=100*kilocalorie_per_mole, 
                          non_aditive_exponent=2.5,
                          qbias_min_seq_sep=3, 
                          qbias_max_seq_sep=np.inf, 
                          qbias_contact_threshold=0.8*nanometers, 
                          forceGroup=4):
    k_qbias = k_qbias.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    qbias = CustomCVForce(f"0.5*{k_qbias}*(q-{q0})^2")
    
    # create bonds
    structure_interactions = read_reference_structure_for_q_calculation_3(oa, 
    reference_pdb_file, reference_chain_name=reference_chain_name,
        min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold, Qflag=0)
    # print(len(structure_interactions))
    # print(structure_interactions)
    # create bond force for q calculation
    normalization = len(structure_interactions)**(1/non_aditive_exponent)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")

    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(forceGroup)
    
    qbias.addCollectiveVariable("q", qvalue)
    # qbias.addGlobalParameter("k_qbias", k_qbias)
    # qbias.addGlobalParameter("q0", q0)
    qbias.setForceGroup(forceGroup)
    return qbias
