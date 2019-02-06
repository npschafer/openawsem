from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
                            'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
                            'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
                            'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}

def read_gamma(gammaFile):
    data = np.loadtxt(gammaFile)
    gamma_direct = data[:210]
    gamma_mediated = data[210:]
    return gamma_direct, gamma_mediated

def inWhichChain(residueId, chain_ends):
    chain_table = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i, end_of_chain_resId in enumerate(chain_ends):
        if end_of_chain_resId < residueId:
            pass
        else:
            return chain_table[i]

def contact_term(oa, k_contact=4.184, z_dependent=False, z_m=1.5, inMembrane=False):
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

    k_relative_mem = 400  # adjust the relative strength of gamma
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
        contact.addComputedValue("alphaMembrane", f"0.5*tanh({eta_switching}*(z+{z_m}))+0.5*tanh({eta_switching}*({z_m}-z))", CustomGBForce.SingleParticle)
        # contact.addComputedValue("alphaMembrane", f"z", CustomGBForce.SingleParticle)
        # contact.addComputedValue("isInMembrane", f"z", CustomGBForce.SingleParticle)
        # contact.addComputedValue("isInMembrane", f"step({z_m}-abs(z))", CustomGBForce.SingleParticle)
        # mediated term
        contact.addEnergyTerm("(1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part;\
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
        contact.addEnergyTerm("(1-alphaMembrane1*alphaMembrane2)*water_part+alphaMembrane1*alphaMembrane2*membrane_part;\
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


def contact_test_term(oa, k_contact=4.184, z_dependent=False, z_m=1.5):
    contact = CustomGBForce()
    gamma_ijm = np.zeros((2, 20, 20))
    contact.addTabulatedFunction("gamma_ijm", Discrete3DFunction(2, 20, 20, gamma_ijm.T.flatten()))
    contact.addComputedValue("rho", f"1", CustomGBForce.ParticlePair)
    contact.addComputedValue("alpha", f"z", CustomGBForce.SingleParticle)
    for i in range(oa.natoms):
        contact.addParticle()
    return contact

'''
# for debug purpose

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

'''
