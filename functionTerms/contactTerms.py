from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
from simtk.unit.quantity import Quantity
import pandas as pd
from Bio.PDB.Polypeptide import three_to_one


gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                            'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                            'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                            'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}


def three_to_index(resName):
    return gamma_se_map_1_letter[three_to_one(resName)]
def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated


def inWhichChain(residueId, chain_ends):
    chain_table = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    for i, end_of_chain_resId in enumerate(chain_ends):
        if end_of_chain_resId < residueId:
            pass
        else:
            return chain_table[i]


def contact_term(oa, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False, membrane_center=0*angstrom, k_relative_mem=1.0, periodic=False, parametersLocation=".", burialPartOn=True, withExclusion=True, forceGroup=22,
                gammaName="gamma.dat", burialGammaName="burial_gamma.dat", membraneGammaName="membrane_gamma.dat", r_min=0.45):
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

    # r_min = .45
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
    gamma_direct, gamma_mediated = read_gamma(os.path.join(parametersLocation, gammaName))

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]
    burial_gamma = np.loadtxt(os.path.join(parametersLocation, burialGammaName))

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
        mem_gamma_direct, mem_gamma_mediated = read_gamma(os.path.join(parametersLocation, membraneGammaName))
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

    contact.addComputedValue("rho", f"isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)))", CustomGBForce.ParticlePair)

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
        contact.addEnergyTerm(f"isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
                                water_part=-res_table(0, resId1, resId2)*{k_contact}*\
                                (gamma_ijm(0, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(0, resName1, resName2)));\
                                membrane_part=-res_table(1, resId1, resId2)*{k_contact}*\
                                (gamma_ijm(1, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(1, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(1, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)));\
                                thetaII=0.25*(1+tanh({eta}*(r-{r_minII})))*(1+tanh({eta}*({r_maxII}-r)));\
                                sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})))",
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
        #                         thetaII=0.25*(1+tanh(eta*(r-{r_minII})))*(1+tanh(eta*({r_maxII}-r)));\
        #                         sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-rho_0)))*(1-tanh({eta_sigma}*(rho2-rho_0)))",
        #                         CustomGBForce.ParticlePair)
        # # direct term
        # contact.addEnergyTerm("isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
        #                         water_part=-res_table(0, resId1, resId2)*k_contact*\
        #                         gamma_ijm(0, resName1, resName2)*theta;\
        #                         membrane_part=-res_table(1, resId1, resId2)*k_contact*\
        #                         gamma_ijm(1, resName1, resName2)*theta;\
        #                         theta=0.25*(1+tanh(eta*(r-r_min)))*(1+tanh(eta*(r_max-r)))",
        #                         CustomGBForce.ParticlePair)
    else:
        # mediated and direct term (write separately may lead to bug)
        contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*{k_contact}*\
                                (gamma_ijm({inMembrane}, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)));\
                                thetaII=0.25*(1+tanh({eta}*(r-{r_minII})))*(1+tanh({eta}*({r_maxII}-r)));\
                                sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})))",
                                CustomGBForce.ParticlePair)
        # # mediated term
        # contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2));\
        #                         sigma_protein=1-sigma_water;\
        #                         thetaII=0.25*(1+tanh(eta*(r-r_minII)))*(1+tanh(eta*(r_maxII-r)));\
        #                         sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})))",
        #                         CustomGBForce.ParticlePair)
        # # direct term
        # contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
        #                         gamma_ijm({inMembrane}, resName1, resName2)*theta;\
        #                         theta=0.25*(1+tanh(eta*(r-r_min)))*(1+tanh(eta*(r_max-r)))",
        #                         CustomGBForce.ParticlePair)

    if burialPartOn:
        # burial term
        for i in range(3):
            contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            contact.addEnergyTerm(f"-0.5*isCb*{k_contact}*burial_gamma_ij(resName, {i})*\
                                        (tanh({burial_kappa}*(rho-rho_min_{i}))+\
                                        tanh({burial_kappa}*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

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

def index_based_contact_term(oa, gamma_folder="ff_contact", k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False, membrane_center=0*angstrom, k_relative_mem=1.0, periodic=False, parametersLocation=".", burialPartOn=True, withExclusion=True, forceGroup=22):
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
    gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    water_gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    protein_gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    # read in seq data.
    seq = oa.seq
    # read in gamma info
    f_direct = np.loadtxt(f"{gamma_folder}/direct.dat")
    f_water = np.loadtxt(f"{gamma_folder}/water.dat")
    f_protein = np.loadtxt(f"{gamma_folder}/protein.dat")
    f_burial = np.loadtxt(f"{gamma_folder}/burial.dat")

    gamma_ijm[0] = f_direct
    water_gamma_ijm[0] = f_water
    protein_gamma_ijm[0] = f_protein
    burial_gamma = f_burial

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]

    k_relative_mem = k_relative_mem  # adjust the relative strength of gamma
    inMembrane = int(inMembrane)
    contact = CustomGBForce()

    m = 0  # water environment
    # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
    res_table = np.zeros((nwell, oa.nres, oa.nres))
    for i in range(oa.nres):
        for j in range(oa.nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                res_table[m][i][j] = 1
            else:
                res_table[m][i][j] = 0


    if z_dependent or inMembrane:
        m = 1  # membrane environment
        f_direct_mem = np.loadtxt(f"{gamma_folder}/direct_mem.dat")
        f_water_mem = np.loadtxt(f"{gamma_folder}/water_mem.dat")
        f_protein_mem = np.loadtxt(f"{gamma_folder}/protein_mem.dat")
        gamma_ijm[m] = f_direct_mem
        water_gamma_ijm[m] = f_water_mem
        protein_gamma_ijm[m] = f_protein_mem
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

    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, oa.nres, oa.nres, gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, oa.nres, oa.nres, water_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, oa.nres, oa.nres, protein_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(oa.nres, 3, burial_gamma.T.flatten()))
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
                                (gamma_ijm(0, resId1, resId2)*theta+thetaII*(sigma_water*water_gamma_ijm(0, resId1, resId2)+\
                                sigma_protein*protein_gamma_ijm(0, resId1, resId2)));\
                                membrane_part=-res_table(1, resId1, resId2)*k_contact*\
                                (gamma_ijm(1, resId1, resId2)*theta+thetaII*(sigma_water*water_gamma_ijm(1, resId1, resId2)+\
                                sigma_protein*protein_gamma_ijm(1, resId1, resId2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)));\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)
    else:
        # mediated and direct term (write separately may lead to bug)
        contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
                                (gamma_ijm({inMembrane}, resId1, resId2)*theta+thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resId1, resId2)+\
                                sigma_protein*protein_gamma_ijm({inMembrane}, resId1, resId2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)));\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)


    if burialPartOn:
        # burial term
        for i in range(3):
            contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            contact.addEnergyTerm(f"-0.5*isCb*k_contact*burial_gamma_ij(resId, {i})*\
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



def expand_contact_table_contact_term(oa, k_contact=4.184, periodic=False, pre=None):
    k_contact *= oa.k_awsem
    # combine direct, burial, mediated.
    # default membrane thickness 1.5 nm

    r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 10  # means j-i > 9
    min_sequence_separation_mem = 13
    # nwell = 16
    eta_switching = 10
    # gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    # water_gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    # protein_gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))

    # read in seq data.
    seq = oa.seq
    # read in gamma info
    if pre is None:
        pre = "expand_contact"
    f_direct = -np.load(f"{pre}/direct.npy")
    f_water = -np.load(f"{pre}/water.npy")
    f_protein = -np.load(f"{pre}/protein.npy")
    f_burial = -np.load(f"{pre}/burial.npy")

    # print("shape of direct", f_direct.shape, len(f_direct.T.flatten()))

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]

    # k_relative_mem = k_relative_mem  # adjust the relative strength of gamma
    inMembrane = 0
    contact = CustomGBForce()

    # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
    res_table = np.zeros((2, oa.nres, oa.nres))
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

    # Discrete3DFunction
    # the tabulated values of the function f(x,y,z), ordered so that values[i+xsize*j+xsize*ysize*k] = f(i,j,k). This must be of length xsize*ysize*zsize.
    contact.addTabulatedFunction("gamma_ijm", Discrete2DFunction(80, 80, f_direct.T.flatten()))
    contact.addTabulatedFunction("water_gamma_ijm", Discrete2DFunction(80, 80, f_water.T.flatten()))
    contact.addTabulatedFunction("protein_gamma_ijm", Discrete2DFunction(80, 80, f_protein.T.flatten()))
    contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(3, 80, f_burial.T.flatten()))
    contact.addTabulatedFunction("res_table", Discrete3DFunction(2, oa.nres, oa.nres, res_table.T.flatten()))

    contact.addPerParticleParameter("resName_with_neighbor")
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

    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    # print(oa.natoms, len(oa.resi), oa.resi, seq)
    for i in range(oa.natoms):
        index = oa.resi[i]
        seq_pre, seq_post = get_pre_and_post(seq, index)
        neighborRes = get_neighbor_res_type(seq_pre, seq_post)
        resName_with_neighbor = int(neighborRes*20 + gamma_se_map_1_letter[seq[oa.resi[i]]])
        contact.addParticle([resName_with_neighbor, gamma_se_map_1_letter[seq[oa.resi[i]]], oa.resi[i], int(i in cb_fixed)])


    # mediated and direct term (write separately may lead to bug)
    contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
                            (gamma_ijm(resName_with_neighbor1, resName_with_neighbor2)*theta+thetaII*(sigma_water*water_gamma_ijm(resName_with_neighbor1, resName_with_neighbor2)+\
                            sigma_protein*protein_gamma_ijm(resName_with_neighbor1, resName_with_neighbor2)));\
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
    # direct term
    # contact.addEnergyTerm(f"-isCb1*isCb2*res_table(0, resId1, resId2)*k_contact*\
    #                         gamma_ijm(0, resId1, resId2)*theta;\
    #                         theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
    #                         CustomGBForce.ParticlePair)

    # burial term
    for i in range(3):
        contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
        contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
    for i in range(3):
        contact.addEnergyTerm(f"-0.5*isCb*k_contact*burial_gamma_ij({i}, resName_with_neighbor)*\
                                    (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                    tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

    print("Number of atom: ", oa.natoms, "Number of residue: ", len(cb_fixed))
    # print(len(none_cb_fixed), len(cb_fixed))
    # for e1 in none_cb_fixed:
    #     for e2 in none_cb_fixed:
    #         if e1 > e2:
    #             continue
    #         contact.addExclusion(e1, e2)
    # for e1 in none_cb_fixed:
    #     for e2 in cb_fixed:
    #         contact.addExclusion(e1, e2)

    # contact.setCutoffDistance(1.1)
    if periodic:
        contact.setNonbondedMethod(contact.CutoffPeriodic)
    else:
        contact.setNonbondedMethod(contact.CutoffNonPeriodic)
    print("Contact cutoff ", contact.getCutoffDistance())
    print("NonbondedMethod: ", contact.getNonbondedMethod())
    contact.setForceGroup(18)
    return contact


def contact_test_term(oa, k_contact=4.184, z_dependent=False, z_m=1.5):
    contact = CustomGBForce()
    gamma_ijm = np.zeros((2, 20, 20))
    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(2, 20, 20, gamma_ijm.T.flatten()))
    contact.addComputedValue("rho", f"1", CustomGBForce.ParticlePair)
    contact.addComputedValue("alpha", f"z", CustomGBForce.SingleParticle)
    for i in range(oa.natoms):
        contact.addParticle()
    return contact

def get_pre_and_post(seq, index):
    n = len(seq)
    if index == 0:
        return seq[0], seq[1]
    elif index == n - 1:
        return seq[index-1], seq[index]
    else:
        return seq[index-1], seq[index+1]

res_type_map_HP = {
    'C': 0,
    'M': 0,
    'F': 0,
    'I': 0,
    'L': 0,
    'V': 0,
    'W': 0,
    'Y': 0,
    'A': 1,
    'H': 1,
    'T': 1,
    'G': 1,
    'P': 1,
    'D': 1,
    'E': 1,
    'N': 1,
    'Q': 1,
    'R': 1,
    'K': 1,
    'S': 1
}

def get_neighbor_res_type(res_pre, res_post):
    table = np.zeros((2,2))
    table[0][0] = 0
    table[0][1] = 1
    table[1][0] = 2
    table[1][1] = 3
    r1 = res_type_map_HP[res_pre]
    r2 = res_type_map_HP[res_pre]
    return int(table[r1][r2])


def hybrid_contact_term(oa, k_contact=4.184, z_m=1.5, membrane_center=0*angstrom, periodic=False, hybrid_gamma_file="hybrid_contact_gamma.dat"):
    k_contact *= oa.k_awsem
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    # combine direct, burial, mediated.
    # default membrane thickness 1.5 nm

    r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 10  # means j-i > 9
    # min_sequence_separation_mem = 10
    nwell = 2
    eta_switching = 10

    # read in seq data.
    seq = oa.seq
    # read in gamma info
    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]

    contact = CustomGBForce()
    gamma = -np.loadtxt(hybrid_gamma_file)
    c = 0
    nwell = 2
    gamma_ijm = np.zeros((nwell, 20, 20))
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))
    burial_gamma_ij = np.zeros((20, 3))
    membrane_burial_gamma_ij = np.zeros((20, 3))
    for ii in range(2):
        for jj in range(3):
            for i in range(20):
                for j in range(i, 20):
                    if jj == 0:
                        gamma_ijm[ii][i][j] = gamma[c]
                        gamma_ijm[ii][j][i] = gamma[c]
                    if jj == 1:
                        protein_gamma_ijm[ii][i][j] = gamma[c]
                        protein_gamma_ijm[ii][j][i] = gamma[c]
                    if jj == 2:
                        water_gamma_ijm[ii][i][j] = gamma[c]
                        water_gamma_ijm[ii][j][i] = gamma[c]
                    c += 1
    for ii in range(2):
        for i in range(3):
            for j in range(20):
                if ii == 0:
                    burial_gamma_ij[j][i] = gamma[c]
                if ii == 1:
                    membrane_burial_gamma_ij[j][i] = gamma[c]
                c += 1

    # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
    for m in range(2):
        res_table = np.zeros((nwell, oa.nres, oa.nres))
        for i in range(oa.nres):
            for j in range(oa.nres):
                resId1 = i
                chain1 = inWhichChain(resId1, oa.chain_ends)
                resId2 = j
                chain2 = inWhichChain(resId2, oa.chain_ends)
                if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                    res_table[m][i][j] = 1
                else:
                    res_table[m][i][j] = 0


    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, 20, 20, gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, 20, 20, water_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, 20, 20, protein_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(20, 3, burial_gamma_ij.T.flatten()))
    contact.addTabulatedFunction("membrane_burial_gamma_ij", Discrete2DFunction(20, 3, membrane_burial_gamma_ij.T.flatten()))
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


    # burial term
    for i in range(3):
        contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
        contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
    for i in range(3):
        contact.addEnergyTerm(f"-0.5*isCb*k_contact*(1 - alphaMembrane)*burial_gamma_ij(resName, {i})*\
                                    (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                    tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)
        contact.addEnergyTerm(f"-0.5*isCb*k_contact*alphaMembrane*membrane_burial_gamma_ij(resName, {i})*\
                                    (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                    tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

    print("Number of atom: ", oa.natoms, "Number of residue: ", len(cb_fixed))
    # print(len(none_cb_fixed), len(cb_fixed))
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
    contact.setForceGroup(18)
    return contact

def disulfide_bond_term(oa, k=1*kilocalorie_per_mole, cutoff=4.2*angstrom, k_bin=100, step_k_bin=20, rho_max=2.2, rho_near=0.2, periodic=False, withExclusion=True, forceGroup=31):
    print("Disulfide Bond term on")
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    cutoff = cutoff.value_in_unit(nanometer)
    # step_k_bin and k_bin in unit of nm^-1.
    k_disulfide_bond = k * oa.k_awsem

    disulfide_bond = CustomGBForce()
    disulfide_bond.addPerParticleParameter("resNameA")
    disulfide_bond.addPerParticleParameter("cysResId")
    disulfide_bond.addPerParticleParameter("cysCB")
    # disulfide_bond.addComputedValue("rho", f"isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)))", CustomGBForce.ParticlePair)
    disulfide_bond.addComputedValue("rhoCys", f"cysCB1*cysCB2*step(abs(cysResId1-cysResId2)-2)*0.5*(1-tanh({k_bin}*(r-{cutoff})))", CustomGBForce.ParticlePair)
    # disulfide_bond.addComputedValue("dummy", f"0", CustomGBForce.SingleParticle)

    cysCount = 0
    for i in range(oa.natoms):
        if oa.seq[oa.resi[i]] == "C" and i in oa.cb:
            isCysCB = 1
        else:
            isCysCB = 0
        cysCount += isCysCB
        disulfide_bond.addParticle([gamma_se_map_1_letter[oa.seq[oa.resi[i]]], oa.resi[i], isCysCB])
    print(f"number of CYS: {cysCount}")

    # disulfide_bond.addEnergyTerm(f"{k_disulfide_bond}*isCb1*isCb2*step(0.5-abs(rho1-rho2))*step(2.5-rho1-rho2)*0.5*(tanh({k_bin}*(r-{cutoff}))-1)",
    #                             CustomGBForce.ParticlePair)
    disulfide_bond.addEnergyTerm(f"{k_disulfide_bond}*cysCB1*cysCB2*min_sep*stepNear*stepSmall*0.5*(tanh({k_bin}*(r-{cutoff}))-1);\
                                    stepNear=0.5*(tanh({step_k_bin}*({rho_near}-abs(rhoCys1-rhoCys2)))+1);\
                                    stepSmall=0.5*(tanh({step_k_bin}*({rho_max}-rhoCys1-rhoCys2))+1);\
                                    min_sep=step(abs(cysResId1-cysResId2)-2)",
                                CustomGBForce.ParticlePair)



    # disulfide_bond.addEnergyTerm(f"{k_disulfide_bond}*isCb1*isCb2*0.5*(tanh({k_bin}*(r-{cutoff}))-1)",
    #                             CustomGBForce.ParticlePair)
    # disulfide_bond.addEnergyTerm(f"{k_disulfide_bond}*isCb1*isCb2", CustomGBForce.ParticlePair)
    # disulfide_bond.addEnergyTerm(f"{k_disulfide_bond}*isCb1*isCb2*(rho1+rho2)", CustomGBForce.ParticlePair)
    if periodic:
        disulfide_bond.setNonbondedMethod(disulfide_bond.CutoffPeriodic)
    else:
        disulfide_bond.setNonbondedMethod(disulfide_bond.CutoffNonPeriodic)
    # disulfide_bond.setCutoffDistance(10)

    # withExclusion won't affect the result. But may speed up the calculation with CPU but slows down for GPU.
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    if withExclusion:
        for e1 in none_cb_fixed:
            for e2 in none_cb_fixed:
                if e1 > e2:
                    continue
                disulfide_bond.addExclusion(e1, e2)
        for e1 in none_cb_fixed:
            for e2 in cb_fixed:
                disulfide_bond.addExclusion(e1, e2)
    print("Disulfide Bond cutoff ", disulfide_bond.getCutoffDistance())
    print("Disulfide NonbondedMethod: ", disulfide_bond.getNonbondedMethod())
    disulfide_bond.setForceGroup(forceGroup)
    return disulfide_bond



def contact_term_shift_well_center(oa, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False, membrane_center=0*angstrom, k_relative_mem=1.0, periodic=False, parametersLocation=".", burialPartOn=True, withExclusion=True, forceGroup=22,
                gammaName="gamma.dat", burialGammaName="burial_gamma.dat", membraneGammaName="membrane_gamma.dat", r_min=0.45, wellCenter=None):
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

    # r_min = .45
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
    gamma_direct, gamma_mediated = read_gamma(os.path.join(parametersLocation, gammaName))

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]
    burial_gamma = np.loadtxt(os.path.join(parametersLocation, burialGammaName))

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
        mem_gamma_direct, mem_gamma_mediated = read_gamma(os.path.join(parametersLocation, membraneGammaName))
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
    if wellCenter:
        # wellCenterInfo = pd.read_csv("/Users/weilu/opt/parameters/side_chain/cbd_cbd_real_contact_symmetric.csv")
        wellCenterInfo = pd.read_csv(wellCenter)
        r_min_direct_table = np.zeros((20, 20))
        r_min_direct_table[three_to_index("GLY"), :] = 2.5
        r_min_direct_table[:, three_to_index("GLY")] = 2.5
        r_max_direct_table = np.zeros((20, 20))
        r_max_direct_table[three_to_index("GLY"), :] = 6.5
        r_max_direct_table[:, three_to_index("GLY")] = 6.5

        r_min_mediated_table = np.zeros((20, 20))
        r_min_mediated_table[three_to_index("GLY"), :] = 6.5
        r_min_mediated_table[:, three_to_index("GLY")] = 6.5
        r_max_mediated_table = np.zeros((20, 20))
        r_max_mediated_table[three_to_index("GLY"), :] = 9.5
        r_max_mediated_table[:, three_to_index("GLY")] = 9.5

        for i, line in wellCenterInfo.iterrows():
            res1 = line["ResName1"]
            res2 = line["ResName2"]
            r_min_ = line["r_min"]
            r_max_ = line["r_max"]
            r_min_direct_table[three_to_index(res1)][three_to_index(res2)] = r_min_ - 0.5
            r_min_direct_table[three_to_index(res2)][three_to_index(res1)] = r_min_ - 0.5

            r_max_direct_table[three_to_index(res1)][three_to_index(res2)] = r_max_ + 1.5
            r_max_direct_table[three_to_index(res2)][three_to_index(res1)] = r_max_ + 1.5

            r_min_mediated_table[three_to_index(res1)][three_to_index(res2)] = r_max_ + 1.5
            r_min_mediated_table[three_to_index(res2)][three_to_index(res1)] = r_max_ + 1.5

            r_max_mediated_table[three_to_index(res1)][three_to_index(res2)] = r_max_ + 4.5
            r_max_mediated_table[three_to_index(res2)][three_to_index(res1)] = r_max_ + 4.5
        # convert unit
        r_min_direct_table /= 10.0
        r_max_direct_table /= 10.0
        r_min_mediated_table /= 10.0
        r_max_mediated_table /= 10.0
        contact.addTabulatedFunction("r_min_direct_table", Discrete2DFunction(20, 20, r_min_direct_table.T.flatten()))
        contact.addTabulatedFunction("r_max_direct_table", Discrete2DFunction(20, 20, r_max_direct_table.T.flatten()))
        contact.addTabulatedFunction("r_min_mediated_table", Discrete2DFunction(20, 20, r_min_mediated_table.T.flatten()))
        contact.addTabulatedFunction("r_max_mediated_table", Discrete2DFunction(20, 20, r_max_mediated_table.T.flatten()))
    contact.addPerParticleParameter("resName")
    contact.addPerParticleParameter("resId")
    contact.addPerParticleParameter("isCb")

    if wellCenter:
        # contact.addComputedValue("rho", f"isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)))", CustomGBForce.ParticlePair)
        contact.addComputedValue("rho", f"isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh({eta}*(r-r_min_direct_table(resName1,resName2))))*(1+tanh({eta}*(r_max_direct_table(resName1,resName2)-r)))", CustomGBForce.ParticlePair)
    else:
        contact.addComputedValue("rho", f"isCb1*isCb2*step(abs(resId1-resId2)-2)*0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)))", CustomGBForce.ParticlePair)

    # if z_dependent:
    #     contact.addComputedValue("isInMembrane", f"step({z_m}-abs(z))", CustomGBForce.SingleParticle)
    # else:
    #     contact.addComputedValue("isInMembrane", "0", CustomGBForce.SingleParticle)


    # contact.addComputedValue("isInMembrane", "1", CustomGBForce.SingleParticle)
    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x, y in zip(oa.cb, oa.ca)]
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
        energy_term = f"isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
                                water_part=-res_table(0, resId1, resId2)*{k_contact}*\
                                (gamma_ijm(0, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(0, resName1, resName2)));\
                                membrane_part=-res_table(1, resId1, resId2)*{k_contact}*\
                                (gamma_ijm(1, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(1, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(1, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)));\
                                thetaII=0.25*(1+tanh({eta}*(r-{r_minII})))*(1+tanh({eta}*({r_maxII}-r)));\
                                sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})))"
        if wellCenter:
            energy_term = f"isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
                                water_part=-res_table(0, resId1, resId2)*{k_contact}*\
                                (gamma_ijm(0, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(0, resName1, resName2)));\
                                membrane_part=-res_table(1, resId1, resId2)*{k_contact}*\
                                (gamma_ijm(1, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm(1, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(1, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh({eta}*(r-r_min)))*(1+tanh({eta}*(r_max-r)));\
                                thetaII=0.25*(1+tanh({eta}*(r-r_minII)))*(1+tanh({eta}*(r_maxII-r)));\
                                sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})));\
                                r_min=r_min_direct_table(resName1,resName2);\
                                r_max=r_max_direct_table(resName1,resName2);\
                                r_minII=r_min_mediated_table(resName1,resName2);\
                                r_maxII=r_max_mediated_table(resName1,resName2);"
        contact.addEnergyTerm(energy_term, CustomGBForce.ParticlePair)
        # # mediated term
        # contact.addEnergyTerm("isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
        #                         water_part=-res_table(0, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm(0, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm(0, resName1, resName2));\
        #                         membrane_part=-res_table(1, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm(1, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm(1, resName1, resName2));\
        #                         sigma_protein=1-sigma_water;\
        #                         thetaII=0.25*(1+tanh(eta*(r-{r_minII})))*(1+tanh(eta*({r_maxII}-r)));\
        #                         sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-rho_0)))*(1-tanh({eta_sigma}*(rho2-rho_0)))",
        #                         CustomGBForce.ParticlePair)
        # # direct term
        # contact.addEnergyTerm("isCb1*isCb2*((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
        #                         water_part=-res_table(0, resId1, resId2)*k_contact*\
        #                         gamma_ijm(0, resName1, resName2)*theta;\
        #                         membrane_part=-res_table(1, resId1, resId2)*k_contact*\
        #                         gamma_ijm(1, resName1, resName2)*theta;\
        #                         theta=0.25*(1+tanh(eta*(r-r_min)))*(1+tanh(eta*(r_max-r)))",
        #                         CustomGBForce.ParticlePair)
    else:
        # mediated and direct term (write separately may lead to bug)
        energy_term = f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*{k_contact}*\
                                (gamma_ijm({inMembrane}, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh({eta}*(r-{r_min})))*(1+tanh({eta}*({r_max}-r)));\
                                thetaII=0.25*(1+tanh({eta}*(r-{r_minII})))*(1+tanh({eta}*({r_maxII}-r)));\
                                sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})))"
        if wellCenter:
            energy_term = f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*{k_contact}*\
                                (gamma_ijm({inMembrane}, resName1, resName2)*theta+thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
                                sigma_protein=1-sigma_water;\
                                theta=0.25*(1+tanh({eta}*(r-r_min)))*(1+tanh({eta}*(r_max-r)));\
                                thetaII=0.25*(1+tanh({eta}*(r-r_minII)))*(1+tanh({eta}*(r_maxII-r)));\
                                sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})));\
                                r_min=r_min_direct_table(resName1,resName2);\
                                r_max=r_max_direct_table(resName1,resName2);\
                                r_minII=r_min_mediated_table(resName1,resName2);\
                                r_maxII=r_max_mediated_table(resName1,resName2);"
            # # direct
            # energy_term = f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*{k_contact}*\
            #                     gamma_ijm({inMembrane}, resName1, resName2)*theta;\
            #                     theta=0.25*(1+tanh({eta}*(r-r_min)))*(1+tanh({eta}*(r_max-r)));\
            #                     r_min=r_min_direct_table(resName1,resName2);\
            #                     r_max=r_max_direct_table(resName1,resName2);"
            # # medaited
            # energy_term = f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*{k_contact}*\
            #                     (thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
            #                     sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
            #                     sigma_protein=1-sigma_water;\
            #                     thetaII=0.25*(1+tanh({eta}*(r-r_minII)))*(1+tanh({eta}*(r_maxII-r)));\
            #                     sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})));\
            #                     r_minII=r_min_mediated_table(resName1,resName2);\
            #                     r_maxII=r_max_mediated_table(resName1,resName2);"
            # # medaited
            # energy_term = f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*{k_contact}*\
            #                     (thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
            #                     sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2)));\
            #                     sigma_protein=1-sigma_water;\
            #                     thetaII=0.25*(1+tanh({eta}*(r-r_minII)))*(1+tanh({eta}*(r_maxII-r)));\
            #                     sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})));\
            #                     r_minII=0.65;\
            #                     r_maxII=0.95;"
        contact.addEnergyTerm(energy_term, CustomGBForce.ParticlePair)
        # # mediated term
        # contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*thetaII*\
        #                         (sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
        #                         sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2));\
        #                         sigma_protein=1-sigma_water;\
        #                         thetaII=0.25*(1+tanh(eta*(r-r_minII)))*(1+tanh(eta*(r_maxII-r)));\
        #                         sigma_water=0.25*(1-tanh({eta_sigma}*(rho1-{rho_0})))*(1-tanh({eta_sigma}*(rho2-{rho_0})))",
        #                         CustomGBForce.ParticlePair)
        # # direct term
        # contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
        #                         gamma_ijm({inMembrane}, resName1, resName2)*theta;\
        #                         theta=0.25*(1+tanh(eta*(r-r_min)))*(1+tanh(eta*(r_max-r)))",
        #                         CustomGBForce.ParticlePair)

    if burialPartOn:
        # burial term
        for i in range(3):
            contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
            contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
        for i in range(3):
            contact.addEnergyTerm(f"-0.5*isCb*{k_contact}*burial_gamma_ij(resName, {i})*\
                                        (tanh({burial_kappa}*(rho-rho_min_{i}))+\
                                        tanh({burial_kappa}*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

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

    contact.setCutoffDistance(1.2)
    if periodic:
        contact.setNonbondedMethod(contact.CutoffPeriodic)
    else:
        contact.setNonbondedMethod(contact.CutoffNonPeriodic)
    print("Contact cutoff ", contact.getCutoffDistance())
    print("NonbondedMethod: ", contact.getNonbondedMethod())
    contact.setForceGroup(forceGroup)
    return contact

'''
# for debug purpose

def contact_term_2(oa, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False, membrane_center=0, k_relative_mem=1.0):
    k_contact *= oa.k_awsem
    # combine direct, burial, mediated.
    # default membrane thickness 1.5 nm

    r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 10  # means j-i > 9
    min_sequence_separation_mem = 13
    nwell = 2
    eta_switching = 10
    gamma_ijm = np.zeros((nwell, 20, 20))
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))

    # read in seq data.
    seq = oa.seq
    # read in gamma info
    gamma_direct, gamma_mediated = read_gamma("gamma.dat")

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]
    burial_gamma = np.loadtxt("burial_gamma.dat")

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
        mem_gamma_direct, mem_gamma_mediated = read_gamma("membrane_gamma.dat")
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

    contact.addComputedValue("rho", "step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)

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
        # mediated term
        contact.addEnergyTerm("((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
                                water_part=-res_table(0, resId1, resId2)*k_contact*thetaII*\
                                (sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(0, resName1, resName2));\
                                membrane_part=-res_table(1, resId1, resId2)*k_contact*thetaII*\
                                (sigma_water*water_gamma_ijm(1, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm(1, resName1, resName2));\
                                sigma_protein=1-sigma_water;\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)
        # direct term
        contact.addEnergyTerm("((1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part);\
                                water_part=-res_table(0, resId1, resId2)*k_contact*\
                                gamma_ijm(0, resName1, resName2)*theta;\
                                membrane_part=-res_table(1, resId1, resId2)*k_contact*\
                                gamma_ijm(1, resName1, resName2)*theta;\
                                theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
                                CustomGBForce.ParticlePair)
    else:
        # mediated term
        contact.addEnergyTerm(f"-res_table({inMembrane}, resId1, resId2)*k_contact*thetaII*\
                                (sigma_water*water_gamma_ijm({inMembrane}, resName1, resName2)+\
                                sigma_protein*protein_gamma_ijm({inMembrane}, resName1, resName2));\
                                sigma_protein=1-sigma_water;\
                                thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                                sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                                CustomGBForce.ParticlePair)
        # direct term
        contact.addEnergyTerm(f"-res_table({inMembrane}, resId1, resId2)*k_contact*\
                                gamma_ijm({inMembrane}, resName1, resName2)*theta;\
                                theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
                                CustomGBForce.ParticlePair)

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
    for e1 in none_cb_fixed:
        for e2 in none_cb_fixed:
            if e1 > e2:
                continue
            contact.addExclusion(e1, e2)
    for e1 in none_cb_fixed:
        for e2 in cb_fixed:
            contact.addExclusion(e1, e2)

    # contact.setCutoffDistance(1.1)
    contact.setNonbondedMethod(CustomGBForce.CutoffNonPeriodic)
    print("Contact cutoff ", contact.getCutoffDistance())
    print("NonbondedMethod: ", contact.getNonbondedMethod())
    contact.setForceGroup(18)
    return contact


def direct_term(oa, k_direct=4.184*1.5):
    k_direct *= oa.k_awsem
    # print(oa.ca, oa.cb)
    # print(oa.bonds)
    # print(oa.nres)  # print 181 for 2xov
    # print(oa.resi)  # print the rsidues index for each atom
    cb = oa.cb
    # gamma = 1
    r_min = .45
    r_max = .65
    eta = 50  # eta actually has unit of nm^-1.
    min_sequence_separation = 10  # means j-i > 9
    nwell = 1
    gamma_ijm = np.zeros((nwell, 20, 20))
    # read in seq data.
    seq = oa.seq
    # read in gamma info
    gamma_direct, gamma_mediated = read_gamma("gamma.dat")

    direct = CustomNonbondedForce(f"-k_direct*gamma_ijm(0, resName1, resName2)*theta; \
    theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r))); \
    eta={eta}")
    # direct = CustomNonbondedForce(f"-k_direct;")
    # direct = CustomNonbondedForce(f"-k_direct*gamma_ijm(0, resName1, resName2);")
    # direct = CustomNonbondedForce(f"-k_direct*gamma_ijm(0, resName1, resName2)*r;")
    direct.addGlobalParameter("k_direct", k_direct)
    direct.addGlobalParameter("rmin", r_min)
    direct.addGlobalParameter("rmax", r_max)



    # add per-particle parameters
    direct.addPerParticleParameter("resName")

    for i in range(oa.natoms):
        direct.addParticle([gamma_se_map_1_letter[seq[oa.resi[i]]]])


    for m in range(nwell):
        count = 0
        for i in range(20):
            for j in range(i, 20):
                gamma_ijm[m][i][j] = gamma_direct[count][0]
                gamma_ijm[m][j][i] = gamma_direct[count][0]
                count += 1

    direct.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, 20, 20, gamma_ijm.flatten()))


    # direct.addInteractionGroup([x for x in cb if x > 0], [x for x in cb if x > 0])
    # direct.addInteractionGroup([x if x > 0 else y for x,y in zip(cb,oa.ca)], [x if x > 0 else y for x,y in zip(cb,oa.ca)])
    # direct.createExclusionsFromBonds(oa.bonds, 11)
    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(cb,oa.ca)]
    # add interaction that are cutoff away
    # don't use this for multi chain simulation.
    for i, x in enumerate(cb_fixed):
        # print(i, x)
        direct.addInteractionGroup([x], cb_fixed[i+min_sequence_separation:])
    # print(cb)

    direct.setForceGroup(16)
    return direct


def burial_term(oa, k_burial=4.184, fastaFile="FastaFileMissing"):
    k_burial *= oa.k_awsem
    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]
    seq = oa.seq
    eta = 50  # eta actually has unit of nm^-1.
    r_min = .45
    r_max = .65
    burial_gamma = np.loadtxt("burial_gamma.dat")

    # return burial
    # if ( lc->chain_no[i]!=lc->chain_no[j] || abs(lc->res_no[j] - lc->res_no[i])>1 )
    burial = CustomGBForce()

    burial_gamma_ij = np.zeros((20, 3))
    burial.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(20, 3, burial_gamma.T.flatten()))

    burial.addPerParticleParameter("resName")
    burial.addPerParticleParameter("resId")
    burial.addPerParticleParameter("isCb")
    burial.addGlobalParameter("k_burial", k_burial)
    burial.addGlobalParameter("eta", eta)
    burial.addGlobalParameter("burial_kappa", burial_kappa)
    burial.addGlobalParameter("rmin", r_min)
    burial.addGlobalParameter("rmax", r_max)
    index = burial.addComputedValue("rho", "step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)
    # print(burial.getComputedValueParameters(index))

    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    for i in range(oa.natoms):
        burial.addParticle([gamma_se_map_1_letter[seq[oa.resi[i]]], oa.resi[i], int(i in cb_fixed)])
    for i in range(3):
        burial.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
        burial.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
    for i in range(3):
        burial.addEnergyTerm(f"-0.5*isCb*k_burial*burial_gamma_ij(resName, {i})*\
                                    (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                    tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

    # burial.addEnergyTerm("-k_burial*rho", CustomGBForce.SingleParticle)
    # burial.addEnergyTerm("-k_burial", CustomGBForce.SingleParticle)


    # print(len(none_cb_fixed), len(cb_fixed))
    for e1 in none_cb_fixed:
        for e2 in none_cb_fixed:
            if e1 > e2:
                continue
            burial.addExclusion(e1, e2)
    for e1 in none_cb_fixed:
        for e2 in cb_fixed:
            burial.addExclusion(e1, e2)

    burial.setForceGroup(17)
    return burial


def mediated_term(oa, k_mediated=4.184*1.5):
    k_mediated *= oa.k_awsem
    # print(oa.nres)  # print 181 for 2xov
    # print(oa.resi)  # print the rsidues index for each atom
    # gamma = 1
    r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 10  # means j-i > 9
    nwell = 1
    water_gamma_ijm = np.zeros((nwell, 20, 20))
    protein_gamma_ijm = np.zeros((nwell, 20, 20))
    # read in seq data.
    seq = oa.seq
    # read in gamma info
    gamma_direct, gamma_mediated = read_gamma("gamma.dat")

    # mediated = CustomNonbondedForce(f"-k_mediated*densityGamma*theta2; \
    # densityGamma=sigmawater_gamma_ijm(0, resName1, resName2); \
    # theta2=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r))); \
    # eta={eta}")
    # mediated = CustomNonbondedForce(f"rho;")

    mediated = CustomGBForce()

    for m in range(nwell):
        count = 0
        for i in range(20):
            for j in range(i, 20):
                water_gamma_ijm[m][i][j] = gamma_mediated[count][1]
                water_gamma_ijm[m][j][i] = gamma_mediated[count][1]
                count += 1

    for m in range(nwell):
        count = 0
        for i in range(20):
            for j in range(i, 20):
                protein_gamma_ijm[m][i][j] = gamma_mediated[count][0]
                protein_gamma_ijm[m][j][i] = gamma_mediated[count][0]
                count += 1
    mediated.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, 20, 20, water_gamma_ijm.flatten()))
    mediated.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, 20, 20, protein_gamma_ijm.flatten()))

    # residue interaction table (step(abs(resId1-resId2)-min_sequence_separation))
    res_table = np.zeros((oa.nres, oa.nres))
    for i in range(oa.nres):
        for j in range(oa.nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            if abs(resId1-resId2)-min_sequence_separation >= 0 or chain1 != chain2:
                res_table[i][j] = 1
            else:
                res_table[i][j] = 0
    mediated.addTabulatedFunction("res_table", Discrete2DFunction(oa.nres, oa.nres, res_table.T.flatten()))
    mediated.addPerParticleParameter("resName")
    mediated.addPerParticleParameter("resId")
    mediated.addPerParticleParameter("isCb")
    mediated.addGlobalParameter("k_mediated", k_mediated)
    mediated.addGlobalParameter("eta", eta)
    mediated.addGlobalParameter("eta_sigma", eta_sigma)
    mediated.addGlobalParameter("rho_0", rho_0)
    mediated.addGlobalParameter("min_sequence_separation", min_sequence_separation)
    mediated.addGlobalParameter("rmin", r_min)
    mediated.addGlobalParameter("rmax", r_max)
    mediated.addGlobalParameter("rminII", r_minII)
    mediated.addGlobalParameter("rmaxII", r_maxII)

    mediated.addComputedValue("rho", "step(abs(resId1-resId2)-2)*0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))", CustomGBForce.ParticlePair)
    # print(burial.getComputedValueParameters(index))

    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    for i in range(oa.natoms):
        mediated.addParticle([gamma_se_map_1_letter[seq[oa.resi[i]]], oa.resi[i], int(i in cb_fixed)])

    mediated.addEnergyTerm("-res_table(resId1, resId2)*k_mediated*thetaII*\
                            (sigma_water*water_gamma_ijm(0, resName1, resName2)+\
                            sigma_protein*protein_gamma_ijm(0, resName1, resName2));\
                            sigma_protein=1-sigma_water;\
                            thetaII=0.25*(1+tanh(eta*(r-rminII)))*(1+tanh(eta*(rmaxII-r)));\
                            sigma_water=0.25*(1-tanh(eta_sigma*(rho1-rho_0)))*(1-tanh(eta_sigma*(rho2-rho_0)))",
                            CustomGBForce.ParticlePair)
    # print(len(none_cb_fixed), len(cb_fixed))
    for e1 in none_cb_fixed:
        for e2 in none_cb_fixed:
            if e1 > e2:
                continue
            mediated.addExclusion(e1, e2)
    for e1 in none_cb_fixed:
        for e2 in cb_fixed:
            mediated.addExclusion(e1, e2)

    mediated.setForceGroup(18)
    return mediated



def index_based_contact_term(oa, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False, membrane_center=0, k_relative_mem=1.0, periodic=False, pre=None):
    z_dependent = False
    inMembrane = False
    k_contact *= oa.k_awsem
    # combine direct, burial, mediated.
    # default membrane thickness 1.5 nm

    r_min = .45
    r_max = .65
    r_minII = .65
    r_maxII = .95
    eta = 50  # eta actually has unit of nm^-1.
    eta_sigma = 7.0
    rho_0 = 2.6
    min_sequence_separation = 10  # means j-i > 9
    min_sequence_separation_mem = 13
    nwell = 2
    eta_switching = 10
    gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    water_gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))
    protein_gamma_ijm = np.zeros((nwell, oa.nres, oa.nres))

    # read in seq data.
    seq = oa.seq
    # read in gamma info
    if pre is None:
        pre = "ff_contact"
    f_direct = np.loadtxt(f"{pre}/direct.dat")
    f_water = np.loadtxt(f"{pre}/water.dat")
    f_protein = np.loadtxt(f"{pre}/protein.dat")
    f_burial = np.loadtxt(f"{pre}/burial.dat")

    gamma_ijm[0] = f_direct
    water_gamma_ijm[0] = f_water
    protein_gamma_ijm[0] = f_protein
    burial_gamma = f_burial

    burial_kappa = 4.0
    burial_ro_min = [0.0, 3.0, 6.0]
    burial_ro_max = [3.0, 6.0, 9.0]

    k_relative_mem = k_relative_mem  # adjust the relative strength of gamma
    inMembrane = int(inMembrane)
    contact = CustomGBForce()

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

    # Discrete3DFunction
    # the tabulated values of the function f(x,y,z), ordered so that values[i+xsize*j+xsize*ysize*k] = f(i,j,k). This must be of length xsize*ysize*zsize.
    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(nwell, oa.nres, oa.nres, gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("water_gamma_ijm", Discrete3DFunction(nwell, oa.nres, oa.nres, water_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("protein_gamma_ijm", Discrete3DFunction(nwell, oa.nres, oa.nres, protein_gamma_ijm.T.flatten()))
    contact.addTabulatedFunction("burial_gamma_ij", Discrete2DFunction(oa.nres, 3, burial_gamma.T.flatten()))
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

    # replace cb with ca for GLY
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    none_cb_fixed = [i for i in range(oa.natoms) if i not in cb_fixed]
    # print(oa.natoms, len(oa.resi), oa.resi, seq)
    for i in range(oa.natoms):
        contact.addParticle([gamma_se_map_1_letter[seq[oa.resi[i]]], oa.resi[i], int(i in cb_fixed)])


    # mediated and direct term (write separately may lead to bug)
    contact.addEnergyTerm(f"-isCb1*isCb2*res_table({inMembrane}, resId1, resId2)*k_contact*\
                            (gamma_ijm({inMembrane}, resId1, resId2)*theta+thetaII*(sigma_water*water_gamma_ijm({inMembrane}, resId1, resId2)+\
                            sigma_protein*protein_gamma_ijm({inMembrane}, resId1, resId2)));\
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
    # direct term
    # contact.addEnergyTerm(f"-isCb1*isCb2*res_table(0, resId1, resId2)*k_contact*\
    #                         gamma_ijm(0, resId1, resId2)*theta;\
    #                         theta=0.25*(1+tanh(eta*(r-rmin)))*(1+tanh(eta*(rmax-r)))",
    #                         CustomGBForce.ParticlePair)

    # burial term
    for i in range(3):
        contact.addGlobalParameter(f"rho_min_{i}", burial_ro_min[i])
        contact.addGlobalParameter(f"rho_max_{i}", burial_ro_max[i])
    for i in range(3):
        contact.addEnergyTerm(f"-0.5*isCb*k_contact*burial_gamma_ij(resId, {i})*\
                                    (tanh(burial_kappa*(rho-rho_min_{i}))+\
                                    tanh(burial_kappa*(rho_max_{i}-rho)))", CustomGBForce.SingleParticle)

    print("Number of atom: ", oa.natoms, "Number of residue: ", len(cb_fixed))
    # print(len(none_cb_fixed), len(cb_fixed))
    # for e1 in none_cb_fixed:
    #     for e2 in none_cb_fixed:
    #         if e1 > e2:
    #             continue
    #         contact.addExclusion(e1, e2)
    # for e1 in none_cb_fixed:
    #     for e2 in cb_fixed:
    #         contact.addExclusion(e1, e2)

    # contact.setCutoffDistance(1.1)
    if periodic:
        contact.setNonbondedMethod(contact.CutoffPeriodic)
    else:
        contact.setNonbondedMethod(contact.CutoffNonPeriodic)
    print("Contact cutoff ", contact.getCutoffDistance())
    print("NonbondedMethod: ", contact.getNonbondedMethod())
    contact.setForceGroup(28)
    return contact

'''
