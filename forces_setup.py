import os
import sys
from simtk.unit import angstrom
from simtk.unit import kilocalorie_per_mole

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
    # print(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *


def set_up_forces(oa, computeQ=False, submode=0, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = [
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        contact_term(oa),
        # contact_term(oa, z_dependent=True, inMembrane=True, membrane_center=membrane_center),
        beta_term_1(oa),
        beta_term_2(oa),
        beta_term_3(oa),
        pap_term_1(oa),
        pap_term_2(oa),
        # membrane_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center),
        # fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        # fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=True),
    ]
    if computeQ:
        forces.append(rg_term(oa))
        forces.append(q_value(oa, "crystal_structure-cleaned.pdb"))
        forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
    return forces
