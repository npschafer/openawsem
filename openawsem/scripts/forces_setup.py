from openawsem.functionTerms import *
from openawsem.helperFunctions.myFunctions import *

try:
    from openmm.unit import angstrom
    from openmm.unit import kilocalorie_per_mole
except ModuleNotFoundError:
    from simtk.unit import angstrom
    from simtk.unit import kilocalorie_per_mole

def set_up_forces(oa, computeQ=False, submode=-1, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = [
        basicTerms.con_term(oa),
        basicTerms.chain_term(oa),
        basicTerms.chi_term(oa),
        basicTerms.excl_term(oa, periodic=False),
        basicTerms.rama_term(oa),
        basicTerms.rama_proline_term(oa),
        basicTerms.rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        contactTerms.contact_term(oa),
        # for membrane protein simulation use contact_term below.
        # contact_term(oa, z_dependent=True, inMembrane=True, membrane_center=membrane_center, k_relative_mem=3),
        hydrogenBondTerms.beta_term_1(oa),
        hydrogenBondTerms.beta_term_2(oa),
        hydrogenBondTerms.beta_term_3(oa),
        hydrogenBondTerms.pap_term_1(oa),
        hydrogenBondTerms.pap_term_2(oa),
        # membraneTerms.membrane_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center),
        # membraneTerms.membrane_preassigned_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center, zimFile="PredictedZim"),
        # templateTerms.er_term(oa),
        # templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),
        debyeHuckelTerms.debye_huckel_term(oa, chargeFile="charge.txt"),
    ]
    if computeQ:
        forces.append(biasTerms.rg_term(oa))
        forces.append(biasTerms.q_value(oa, "crystal_structure-cleaned.pdb", forceGroup=1))
        # forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", residueIndexGroup=list(range(0, 15)), forceGroup=1))
    if submode == 0:
        additional_forces = [
            # contact_term(oa),
        ]
        forces += additional_forces
    return forces
