import os
import sys
try:
    from openmm.unit import angstrom
    from openmm.unit import kilocalorie_per_mole
except ModuleNotFoundError:
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

# k_rg = 20, 3TM
# g_all = [
#     [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28],
#     [73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93],
#     [136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157],
# ]
# k_rg = 20, 4TM
# g_all = [
#     [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25],
#     [67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88],
#     [89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106],
#     [133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157],
# ]


def cbd_excl_term(oa, k=1*kilocalorie_per_mole, periodic=False, r_excl=0.7, fileLocation='../../cbd_cbd_real_contact_symmetric.csv', forceGroup=24):
    # Cb domain Excluded volume
    # With residue specific parameters
    # a harmonic well with minimum at the database histogram peak.
    # and 1 kT(0.593 kcal/mol) penalty when the distance is at r_min of the database.
    # multiply interaction strength by overall scaling
    # Openawsem doesn't have the distance range (r_excl) change from 0.35 to 0.45 when the sequence separtation more than 5
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_excl = k * oa.k_awsem
    excl = CustomNonbondedForce(f"{k_excl}*step(r_max(res1,res2)-r)*((r-r_max(res1,res2))/(r_max(res1,res2)-r_min(res1,res2)))^2")
    excl.addPerParticleParameter("res")

    gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                                'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                                'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                                'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
    for i in range(oa.natoms):
        excl.addParticle([gamma_se_map_1_letter[oa.seq[oa.resi[i]]]])
    # print(oa.ca)
    # print(oa.bonds)
    # print(oa.cb)

    r_min_table = np.zeros((20,20))
    r_max_table = np.zeros((20,20))
    # fileLocation = '/Users/weilu/Research/server/mar_2020/cmd_cmd_exclude_volume/cbd_cbd_real_contact_symmetric.csv'
    df = pd.read_csv(fileLocation)
    for i, line in df.iterrows():
        res1 = line["ResName1"]
        res2 = line["ResName2"]
        r_min_table[gamma_se_map_1_letter[three_to_one(res1)]][gamma_se_map_1_letter[three_to_one(res2)]] = line["r_min"] / 10.0   # A to nm
        r_min_table[gamma_se_map_1_letter[three_to_one(res2)]][gamma_se_map_1_letter[three_to_one(res1)]] = line["r_min"] / 10.0
        r_max_table[gamma_se_map_1_letter[three_to_one(res1)]][gamma_se_map_1_letter[three_to_one(res2)]] = line["r_max"] / 10.0
        r_max_table[gamma_se_map_1_letter[three_to_one(res2)]][gamma_se_map_1_letter[three_to_one(res1)]] = line["r_max"] / 10.0

    excl.addTabulatedFunction("r_min", Discrete2DFunction(20, 20, r_min_table.T.flatten()))
    excl.addTabulatedFunction("r_max", Discrete2DFunction(20, 20, r_max_table.T.flatten()))
    excl.addInteractionGroup([x for x in oa.cb if x > 0], [x for x in oa.cb if x > 0])

    excl.setCutoffDistance(r_excl)
    if periodic:
        excl.setNonbondedMethod(excl.CutoffPeriodic)
    else:
        excl.setNonbondedMethod(excl.CutoffNonPeriodic)

    # excl.setNonbondedMethod(excl.CutoffNonPeriodic)
    excl.createExclusionsFromBonds(oa.bonds, 1)
    excl.setForceGroup(forceGroup)
    print("cb domain exlcude volume term On")
    return excl

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

def SideToZ_m(side):
    side = side.strip()
    if side == "down":
        return -1.6
    if side == "up":
        return 1.6
    if side == "middle":
        return 0

def membrane_preassigned_side_term(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, zimFile="PredictedZimSide_4TH", forceGroup=24):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_membrane = k * oa.k_awsem
    membrane = CustomExternalForce(f"{k_membrane}*(abs(z-{membrane_center}-z_m))")
    membrane.addPerParticleParameter("z_m")

    with open(zimFile) as f:
        a = f.readlines()

    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    # print(cb_fixed)
    for i in cb_fixed:
        z_m = SideToZ_m(a[oa.resi[i]])
        # print(oa.resi[i])
        membrane.addParticle(i, [z_m])
    membrane.setForceGroup(forceGroup)
    return membrane

def positive_inside_rule_term(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, forceGroup=25):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_positive_inside_rule = k * oa.k_awsem
    positive_inside_rule = CustomExternalForce(f"{k_positive_inside_rule}*((z-{membrane_center}-z_m)^2)")
    positive_inside_rule.addPerParticleParameter("z_m")

    positive_inside_residue_table = {"G":0, "A":0, "V":0, "C":0, "P":0, "L":0, "I":0, "M":0, "W":0, "F":0,
                                    "S":0, "T":0, "Y":0, "N":0, "Q":0,
                                    "K":-1, "R":-1, "H":0,
                                    "D":0, "E":0}
    thickness = 1.5
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    # print(cb_fixed)
    for i in cb_fixed:
        z_m = positive_inside_residue_table[oa.seq[oa.resi[i]]] * thickness
        # print(oa.resi[i])
        positive_inside_rule.addParticle(i, [z_m])
    positive_inside_rule.setForceGroup(forceGroup)
    return positive_inside_rule

def positive_inside_rule_term_v2(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, forceGroup=25):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_positive_inside_rule = k * oa.k_awsem
    positive_inside_rule = CustomExternalForce(f"{k_positive_inside_rule}*step(z-{membrane_center}-z_m)*((z-{membrane_center}-z_m)^2)")
    positive_inside_rule.addPerParticleParameter("z_m")

    positive_inside_residue_table = {"G":0, "A":0, "V":0, "C":0, "P":0, "L":0, "I":0, "M":0, "W":0, "F":0,
                                    "S":0, "T":0, "Y":0, "N":0, "Q":0,
                                    "K":-1, "R":-1, "H":0,
                                    "D":0, "E":0}
    thickness = 1.5
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    # print(cb_fixed)
    for i in cb_fixed:
        z_m = positive_inside_residue_table[oa.seq[oa.resi[i]]] * thickness
        # print(oa.resi[i])
        positive_inside_rule.addParticle(i, [z_m])
    positive_inside_rule.setForceGroup(forceGroup)
    return positive_inside_rule

def positive_inside_rule_term_v3(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, forceGroup=25):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_positive_inside_rule = k * oa.k_awsem
    positive_inside_rule = CustomExternalForce(f"{k_positive_inside_rule}*step(z-{membrane_center}-z_m)*((z-{membrane_center}-z_m)^2)")
    positive_inside_rule.addPerParticleParameter("z_m")

    positive_inside_residue_table = {"G":0, "A":0, "V":0, "C":0, "P":0, "L":0, "I":0, "M":0, "W":0, "F":0,
                                    "S":0, "T":0, "Y":0, "N":0, "Q":0,
                                    "K":-1, "R":-1, "H":0,
                                    "D":0, "E":0}
    thickness = 1.5
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    # print(cb_fixed)
    for i in cb_fixed:
        z_m = positive_inside_residue_table[oa.seq[oa.resi[i]]] * thickness
        # print(oa.resi[i])
        if z_m != 0:
            positive_inside_rule.addParticle(i, [z_m])
    positive_inside_rule.setForceGroup(forceGroup)
    return positive_inside_rule

def set_topology(oa, zimFile="PredictedZimSide"):
    forces = [
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        membrane_preassigned_term(oa, k=10*kilocalorie_per_mole, membrane_center=0*angstrom, zimFile=zimFile),
        # membrane_preassigned_side_term(oa, k=100*kilocalorie_per_mole, membrane_center=0*angstrom, zimFile=zimFile),
        # pulling_term(oa, k_pulling=5*4.184, forceDirect="x", appliedToResidue=7),
        # pulling_term(oa, k_pulling=5*4.184, forceDirect="-x", appliedToResidue=397)
    ]
    return forces


def pull_part(oa, zimFile="PredictedZim"):
    forces = [
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        membrane_preassigned_term(oa, k=10*kilocalorie_per_mole, membrane_center=0*angstrom, zimFile=zimFile),
        pulling_term(oa, k_pulling=5*4.184, forceDirect="x", appliedToResidue="FIRST"),
        pulling_term(oa, k_pulling=5*4.184, forceDirect="-x", appliedToResidue="LAST"),
        # positive_inside_rule_term_v3(oa, k=10*kilocalorie_per_mole, membrane_center=0*angstrom, forceGroup=8)
    ]
    return forces

# residues_has_crystal_no_filled = list(range(36)) + list(range(48, 109)) + list(range(124,158))

def set_up_forces(oa, computeQ=False, submode=0, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = []
    forces = [
        con_term(oa),
        chain_term(oa),
        chi_term(oa),
        excl_term(oa, periodic=False),
        rama_term(oa),
        rama_proline_term(oa),
        rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
        # contact_term(oa),
        contact_term(oa, inMembrane=True),
        # for membrane protein simulation use contact_term below.
        # contact_term(oa, z_dependent=True, inMembrane=True, membrane_center=membrane_center, k_relative_mem=3),
        # beta_term_1(oa),
        # beta_term_2(oa),
        # beta_term_3(oa),
        # pap_term_1(oa),
        # pap_term_2(oa),
        # er_term(oa),
        membrane_term(oa, k=2*kilocalorie_per_mole, membrane_center=membrane_center),
        # membrane_preassigned_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center, zimFile="PredictedZim"),
        # fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=True),
        # debye_huckel_term(oa, chargeFile="charge.txt"),
        # debye_huckel_term(oa)
    ]
    if computeQ:
        forces.append(rg_term(oa))
        forces.append(q_value(oa, "crystal_structure-cleaned.pdb"))
        forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", a=0.2, forceGroup=3))
        forces.append(positive_inside_rule_term_v3(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, forceGroup=8))
        # forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=GlobularPart, forceGroup=4))
        # forces.append(partial_q_value(oa, "crystal_structure-cleaned.pdb", a=0.2, startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=MembranePart, forceGroup=5))

    if submode == 1:
        additional_forces = [
            con_term(oa),
            chain_term(oa),
            chi_term(oa),
            rama_term(oa),
            rama_proline_term(oa),
            rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
            excl_term(oa, periodic=False),
            fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
            helical_term(oa, k_helical=2*4.184, inMembrane=True),
            # positive_inside_rule_term_v2(oa),
            contact_term(oa, k_contact=4.184, z_dependent=True, inMembrane=True, membrane_center=membrane_center,
                                    k_relative_mem=2, periodic=False, parametersLocation=contactParameterLocation, burialPartOn=True),
            membrane_term(oa, k=2*kilocalorie_per_mole, membrane_center=membrane_center),
        ]
        forces += additional_forces
        for i in range(0, oa.nres-8, 2):
            forces.append(single_helix_orientation_bias_term(oa, k=5*kilocalorie_per_mole, z_m=1.5, membrane_center=membrane_center, atomGroup=[i, i+2, i+4, i+6, i+8]))

    if submode == 2:
        # helix rg k = 2
        additional_forces = [
            con_term(oa),
            chain_term(oa),
            chi_term(oa),
            rama_term(oa),
            rama_proline_term(oa),
            rama_ssweight_term(oa, k_rama_ssweight=2*8.368),
            excl_term(oa, periodic=False),
            fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
            helical_term(oa, k_helical=2*4.184, inMembrane=True),
            # positive_inside_rule_term_v2(oa),
            contact_term(oa, k_contact=4.184, z_dependent=True, inMembrane=True, membrane_center=membrane_center,
                                    k_relative_mem=2, periodic=False, parametersLocation=contactParameterLocation, burialPartOn=True),
            membrane_term(oa, k=2*kilocalorie_per_mole, membrane_center=membrane_center),
        ]
        forces += additional_forces
        for i in range(0, oa.nres-8, 2):
            forces.append(single_helix_orientation_bias_term(oa, k=2*kilocalorie_per_mole, z_m=1.5, membrane_center=membrane_center, atomGroup=[i, i+2, i+4, i+6, i+8]))


    return forces
