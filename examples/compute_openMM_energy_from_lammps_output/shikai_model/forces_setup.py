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



def contact_term_check(oa, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False, membrane_center=0*angstrom, k_relative_mem=1.0, periodic=False, parametersLocation=".", burialPartOn=True, withExclusion=True, forceGroup=22):
    if isinstance(k_contact, float) or isinstance(k_contact, int):
        k_contact = k_contact * oa.k_awsem   # just for backward comptable
    elif isinstance(k_contact, Quantity):
        k_contact = k_contact.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
        k_contact = k_contact * oa.k_awsem
    else:
        print(f"Unknown input, {k_contact}, {type(k_contact)}")
    # combine direct, burial, mediated.
    # default membrane thickness 1.5 nm
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm

    r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 10  # means j-i > 9
    min_sequence_separation_mem = 10
    nwell = 2
    eta_switching = 10
    gamma_ijm = np.zeros((nwell, 20, 20))
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))

    # read in seq data.
    seq = oa.seq
    # read in gamma info
    gamma_direct, gamma_mediated = read_gamma(os.path.join(parametersLocation, "gamma.dat"))

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]
    burial_gamma = np.loadtxt(os.path.join(parametersLocation, "burial_gamma.dat"))

    k_relative_mem = k_relative_mem  # adjust the relative strength of gamma
    inMembrane = int(inMembrane)
    contact = CustomGBForce()

    m = 0  # water environment
    count = 0
    for i in range(20):
        for j in range(i, 20):
            gamma_ijm[m][i][j] = gamma_direct[count][0]
            gamma_ijm[m][j][i] = gamma_direct[count][0]
            count += 1
    count = 0
    for i in range(20):
        for j in range(i, 20):
            water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
            water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
            count += 1
    count = 0
    for i in range(20):
        for j in range(i, 20):
            protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
            protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
            count += 1
    # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
    res_table = np.zeros((nwell, oa.nres, oa.nres))
    for i in range(oa.nres):
        for j in range(oa.nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                res_table[0][i][j] = 1
            else:
                res_table[0][i][j] = 0


    if z_dependent or inMembrane:
        mem_gamma_direct, mem_gamma_mediated = read_gamma(os.path.join(parametersLocation, "membrane_gamma.dat"))
        m = 1  # membrane environment
        count = 0
        for i in range(20):
            for j in range(i, 20):
                gamma_ijm[m][i][j] = mem_gamma_direct[count][0]*k_relative_mem
                gamma_ijm[m][j][i] = mem_gamma_direct[count][0]*k_relative_mem
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                water_gamma_ijm[m][i][j] = mem_gamma_mediated[count][1]*k_relative_mem
                water_gamma_ijm[m][j][i] = mem_gamma_mediated[count][1]*k_relative_mem
                count += 1
        count = 0
        for i in range(20):
            for j in range(i, 20):
                protein_gamma_ijm[m][i][j] = mem_gamma_mediated[count][0]*k_relative_mem
                protein_gamma_ijm[m][j][i] = mem_gamma_mediated[count][0]*k_relative_mem
                count += 1
        for i in range(oa.nres):
            for j in range(oa.nres):
                resId1 = i
                chain1 = inWhichChain(resId1, oa.chain_ends)
                resId2 = j
                chain2 = inWhichChain(resId2, oa.chain_ends)
                if abs(resId1-resId2)-min_sequence_separation_mem >= 0 or chain1 != chain2:
                    res_table[m][i][j] = 1
                else:
                    res_table[m][i][j] = 0

    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, 20, 20, gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, 20, 20, water_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, 20, 20, protein_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(20, 3, burial_gamma.T.flatten()))
    contact.addTabulatedFunction("res_table", Discrete3DFunction(nwell, oa.nres, oa.nres, res_table.T.flatten()))

    contact.addPerParticleParameter("resName")
    contact.addPerParticleParameter("resId")
    contact.addPerParticleParameter("isCb")
    contact.addGlobalParameter("k_contact", k_contact)
    contact.addGlobalParameter("eta", eta)
    contact.addGlobalParameter("eta_sigma", eta_sigma)
    contact.addGlobalParameter("rho_0", rho_0)
    contact.addGlobalParameter("min_sequence_separation", min_sequence_separation)
    contact.addGlobalParameter("rmin", r_min)
    contact.addGlobalParameter("rmax", r_max)
    contact.addGlobalParameter("rminII", r_minII)
    contact.addGlobalParameter("rmaxII", r_maxII)
    contact.addGlobalParameter("burial_kappa", burial_kappa)

    contact.addComputedValue("rho", "isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)

    # if z_dependent:
    #     contact.addComputedValue("isInMembrane", f"step({z_m}-abs(z))", CustomGBForce.SingleParticle)
    # else:
    #     contact.addComputedValue("isInMembrane", "0", CustomGBForce.SingleParticle)


    # contact.addComputedValue("isInMembrane", "1", CustomGBForce.SingleParticle)
    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    # print(oa.natoms, len(oa.resi), oa.resi, seq)
    for i in range(oa.natoms):
        contact.addParticle([gamma_se_map_1_letter[seq[oa.resi[i]]], oa.resi[i], int(i in cb_fixed)])


    if z_dependent:
        # print(f"0.5*tanh({eta_switching}*(z+{z_m}))+0.5*tanh({eta_switching}*({z_m}-z))")
        contact.addComputedValue("alphaMembrane", f"0.5*tanh({eta_switching}*((z-{membrane_center})+{z_m}))+0.5*tanh({eta_switching}*({z_m}-(z-{membrane_center})))", CustomGBForce.SingleParticle)
        # contact.addComputedValue("alphaMembrane", f"z", CustomGBForce.SingleParticle)
        # contact.addComputedValue("isInMembrane", f"z", CustomGBForce.SingleParticle)
        # contact.addComputedValue("isInMembrane", f"step({z_m}-abs(z))", CustomGBForce.SingleParticle)

        # mediated and direct term (write separately may lead to bug)
        contact.addEnergyTerm("isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
                                water_part=-res_table(0, resId1, resId2)*k_contact*\
                                (gamma_ijm(0, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(0, resName1, resName2)));\
                                membrane_part=-res_table(1, resId1, resId2)*k_contact*\
                                (gamma_ijm(1, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(1, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(1, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)));\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)
        # # mediated term
        # contact.addEnergyTerm("isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
        #                         water_part=-res_table(0, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm(0, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm(0, resName1, resName2));\
        #                         membrane_part=-res_table(1, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm(1, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm(1, resName1, resName2));\
        #                         sigma_protein=1-sigma_water;\
        #                         thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
        #                         sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
        #                         CustomGBForce.ParticlePair)
        # # direct term
        # contact.addEnergyTerm("isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
        #                         water_part=-res_table(0, resId1, resId2)*k_contact*\
        #                         gamma_ijm(0, resName1, resName2)*theta;\
        #                         membrane_part=-res_table(1, resId1, resId2)*k_contact*\
        #                         gamma_ijm(1, resName1, resName2)*theta;\
        #                         theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
        #                         CustomGBForce.ParticlePair)
    else:
        # mediated and direct term (write separately may lead to bug)
        contact.addEnergyTerm(f"-0.75*isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
                                (gamma_ijm({inMembrane}, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)));\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)
        # # mediated term
        # contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2));\
        #                         sigma_protein=1-sigma_water;\
        #                         thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
        #                         sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
        #                         CustomGBForce.ParticlePair)
        # # direct term
        # contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
        #                         gamma_ijm({inMembrane}, resName1, resName2)*theta;\
        #                         theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
        #                         CustomGBForce.ParticlePair)

    if burialPartOn:
        # burial term
        for i in range(3):
            contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            contact.addEnergyTerm(f"-0.5*isCb*k_contact*burial_gamma_ij(resName, {i})*\
                                        (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                        tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

    print("Number of atom: ", oa.natoms, "Number of residue: ", len(cb_fixed))
    # print(len(none_cb_fixed), len(cb_fixed))

    # withExclusion won't affect the result. But may speed up the calculation with CPU but slows down for GPU.
    if withExclusion:
        for e1 in none_cb_fixed:
            for e2 in none_cb_fixed:
                if e1 > e2:
                    continue
                contact.addExclusion(e1, e2)
        for e1 in none_cb_fixed:
            for e2 in cb_fixed:
                contact.addExclusion(e1, e2)

    # contact.setCutoffDistance(1.1)
    if periodic:
        contact.setNonbondedMethod(contact.CutoffPeriodic)
    else:
        contact.setNonbondedMethod(contact.CutoffNonPeriodic)
    print("Contact cutoff ", contact.getCutoffDistance())
    print("NonbondedMethod: ", contact.getNonbondedMethod())
    contact.setForceGroup(forceGroup)
    return contact

def excl_term_check(oa, k_excl=8368, r_excl=0.35, periodic=False, forceGroup=20):
    # add excluded volume
    # Still need to add element specific parameters
    # 8368 = 20 * 4.184 * 100 kJ/nm^2, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_excl *= oa.k_awsem
    excl = CustomNonbondedForce("k_excl*step(r0-r)*(r-r0)^2")
    excl.addGlobalParameter("k_excl", k_excl)
    excl.addGlobalParameter("r0", r_excl)
    for i in range(oa.natoms):
        excl.addParticle()
    # print(oa.ca)
    # print(oa.bonds)
    # print(oa.cb)
    excl.addInteractionGroup(oa.ca, oa.ca)
    excl.addInteractionGroup([x for x in oa.cb if x > 0], [x for x in oa.cb if x > 0])
    excl.addInteractionGroup(oa.ca, [x for x in oa.cb if x > 0])
    excl.addInteractionGroup(oa.o, oa.o)

    excl.setCutoffDistance(r_excl)
    if periodic:
        excl.setNonbondedMethod(excl.CutoffPeriodic)
    else:
        excl.setNonbondedMethod(excl.CutoffNonPeriodic)

    # excl.setNonbondedMethod(excl.CutoffNonPeriodic)
    excl.createExclusionsFromBonds(oa.bonds, 1)
    excl.setForceGroup(forceGroup)
    return excl

def set_up_forces(oa, computeQ=False, submode=0, contactParameterLocation=".", membrane_center=-0*angstrom):
    # apply forces
    forces = [
        con_term(oa, forceGroup=11, k_con=167360, bond_lengths=[.377, .241, .250, .154]),   # HB used in Lammps version is not included here. To match the energy exactly, you need set the k_coeff of bond 5 to be 0.
        chain_term(oa, forceGroup=12),
        chi_term(oa, forceGroup=13),
        excl_term_check(oa, periodic=False, forceGroup=14, r_excl=0.35),
        rama_term(oa, forceGroup=15, phi_i=[-1.74, -1.265, 1.3], psi_i=[2.138, -0.318, 0.5]),
        rama_proline_term(oa, forceGroup=15),
        rama_ssweight_term(oa, k_rama_ssweight=8.368, forceGroup=15, sigma=[99.0, 15.398], phi_i=[-1.1, -2.25]),
        tbm_q_term(oa, 2*200*4.184)
        # contact_term_check(oa),
        # contact_term(oa, z_dependent=True, inMembrane=True, membrane_center=membrane_center),
        # beta_term_1(oa),
        # beta_term_2(oa),
        # beta_term_3(oa),
        # pap_term_1(oa),
        # pap_term_2(oa),
        # er_term(oa),
        # membrane_term(oa, k=1*kilocalorie_per_mole, membrane_center=membrane_center),
        # fragment_memory_term(oa, frag_file_list_file="./frags.mem.single", npy_frag_table="./frags.npy", UseSavedFragTable=True, min_seq_sep=2, max_seq_sep=9),
        # fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=True),
    ]
    if computeQ:
        forces.append(rg_term(oa))
        forces.append(q_value(oa, "crystal_structure-cleaned.pdb", forceGroup=1))
        forces.append(qc_value(oa, "crystal_structure-cleaned.pdb"))
    return forces
