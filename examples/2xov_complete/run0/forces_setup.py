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

k_rg = 20
g_all = [
    [96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113],
    [139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161],
    [171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186],
    [196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212],
    [223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241],
    [252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266],
]

MembranePart = [91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270]
# GlobularPart = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85]
GlobularPart = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62]


def group_pulling_term(oa, k=1*kilocalorie_per_mole, forceDirect="z", appliedToResidues=-1, membrane_center=0*angstrom):
    # appliedToResidues can be a list of residue index.
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_pulling = k * oa.k_awsem
    pulling = CustomExternalForce(f"{forceDirect}")
    for i in range(oa.natoms):
        if appliedToResidues == -1:
            pulling.addParticle(i)
            n = len(oa.resi)
        elif oa.resi[i] in appliedToResidues:
            pulling.addParticle(i)
            n = len(appliedToResidues)
        # if oa.resi[i] == appliedToResidue:
        #     pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])

    harmonic = CustomCVForce(f"{k_pulling}*(average/{n}-{membrane_center})^2")
    harmonic.addCollectiveVariable("average", pulling)
    harmonic.setForceGroup(29)
    return harmonic

def group_pulling_linear_term(oa, k=1*kilocalorie_per_mole, forceDirect="z", appliedToResidues=-1, membrane_center=0*angstrom):
    # appliedToResidues can be a list of residue index.
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_pulling = k * oa.k_awsem

    pulling = CustomExternalForce(f"{forceDirect}")
    for i in range(oa.natoms):
        if appliedToResidues == -1:
            pulling.addParticle(i)
            n = oa.natoms
        elif oa.resi[i] in appliedToResidues:
            pulling.addParticle(i)
            n = len(appliedToResidues)
        # if oa.resi[i] == appliedToResidue:
        #     pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])

    linear = CustomCVForce(f"{k_pulling}*abs(average/{n} - {membrane_center})")
    linear.addCollectiveVariable("average", pulling)
    linear.setForceGroup(29)
    return linear

def pull_part(oa):
    forces = [
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        membrane_preassigned_term(oa, k=10*kilocalorie_per_mole, membrane_center=0*angstrom, zimFile="PredictedZim"),
        pulling_term(oa, k_pulling=5*4.184, forceDirect="x", appliedToResidue=96),
        pulling_term(oa, k_pulling=5*4.184, forceDirect="-x", appliedToResidue=266)
    ]
    return forces

def set_up_forces(oa, params, computeQ=False, submode=0, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = [
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        # contact_term(oa, k_contact=params.k_contact, z_dependent=params.z_dependent, inMembrane=params.inMembrane,
        #                 k_relative_mem=params.k_relative_mem, periodic=params.periodic),
        # index_based_contact_term(oa, pre="ff_contact/"),
        # expand_contact_table_contact_term(oa, pre="/Users/weilu/Research/server/may_2019/openMM_multiLetter/symmetric"),
        # expand_contact_table_contact_term(oa, pre="../../symmetric"),

        beta_term_1(oa),
        beta_term_2(oa),
        beta_term_3(oa),
        pap_term_1(oa),
        pap_term_2(oa),

        # group_pulling_term(oa, k=10*kilocalorie_per_mole, forceDirect="z", appliedToResidues=-1, membrane_center=membrane_center),
        # group_pulling_linear_term(oa, k=100*kilocalorie_per_mole, forceDirect="z", appliedToResidues=-1, membrane_center=membrane_center),
        # fragment_memory_term(oa, frag_file_list_file="./HA_combined.mem", npy_frag_table="./HA_combined.npy", UseSavedFragTable=True),
        fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=True),

        # fragment_memory_term(oa, frag_file_list_file="./combined.mem", npy_frag_table="./combined.npy", UseSavedFragTable=True),
        # membrane_term(oa, k=2*kilocalorie_per_mole, membrane_center=membrane_center),
        membrane_preassigned_term(oa, k=5*kilocalorie_per_mole, membrane_center=membrane_center, zimFile="PredictedZim"),
        z_dependent_helical_term(oa, k_helical=1*4.184, membrane_center=membrane_center),
        # fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        # fragment_memory_term(oa, frag_file_list_file="./new_frags.mem", npy_frag_table="./new_frags.npy", UseSavedFragTable=True),
        er_term(oa),
        # tbm_q_term(oa, k_tbm_q=2000),
        # helical_term(oa, k_helical=2*4.184, inMembrane=True),
        # membrane_preassigned_term(oa, k_membrane=10*4.184, membrane_center=0),
        # pulling_term(oa, k_pulling=10*4.184, forceDirect="x", appliedToResidue=0),
        # pulling_term(oa, k_pulling=10*4.184, forceDirect="-x", appliedToResidue=180),
        # rg_bias_term(oa, k_rg=100*4.184, rg0=0),
        # single_helix_orientation_bias_term(oa, k_rg=4.184, rg0=0, atomGroup=list(range(11, 46))),
    ]
    for i in range(len(g_all)):
        forces.append(single_helix_orientation_bias_term(oa, k=k_rg*kilocalorie_per_mole, z_m=1.5, membrane_center=membrane_center, atomGroup=g_all[i]))
    if computeQ:
        forces.append(rg_term(oa))
        forces.append(q_value(oa, "crystal_structure-cleaned.pdb"))
        forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
        forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=GlobularPart, forceGroup=4))
        forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", a=0.2, startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=MembranePart, forceGroup=5))
    if submode == 0:
        forces.append(expand_contact_table_contact_term(oa, pre="../../symmetric"))
    elif submode == 1:
        # forces.append(contact_term(oa, k_contact=params.k_contact, z_dependent=False, inMembrane=False,
        #                              k_relative_mem=1, periodic=False))
        forces.append(contact_term(oa, k_contact=params.k_contact, z_dependent=False, inMembrane=True,
                                    k_relative_mem=1, periodic=False, parametersLocation=contactParameterLocation, burialPartOn=False))
    elif submode == 2:
        forces.append(hybrid_contact_term(oa, periodic=False, membrane_center=membrane_center, hybrid_gamma_file="../../hybrid_gamma.dat"))
    elif submode == 3:
        forces.append(contact_term(oa, k_contact=1 * 4.184, z_dependent=True, inMembrane=True, membrane_center=membrane_center,
                                    k_relative_mem=1, periodic=False, parametersLocation=contactParameterLocation, burialPartOn=True))

    return forces
