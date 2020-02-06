from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

se_map_1_letter = {'A': 0,  'P': 1,  'K': 2,  'N': 3,  'R': 4,
                   'F': 5,  'D': 6,  'Q': 7,  'E': 8,  'G': 9,
                   'I': 10, 'H': 11, 'L': 12, 'C': 13, 'M': 14,
                   'S': 15, 'T': 16, 'Y': 17, 'V': 18, 'W': 19}

def isChainStart(residueId, chain_starts, n=2):
    # return true if residue is near chain starts.
    # n=0 means always return False
    # n=1 means only return True if residue is the the first residue of a chain.
    # n=2 means return True if residue is the first or the one nearest to the first residue of a chain.
    atBegin = False
    for i in range(n):
        if (residueId-i) in chain_starts:
            atBegin = True
    return atBegin

def isChainEnd(residueId, chain_ends, n=2):
    # return true if residue is near chain ends.
    # n=0 means always return False
    # n=1 means only return True if residue is the the last residue of a chain.
    # n=2 means return True if residue is the last or the one nearest to the last residue of a chain.
    atEnd = False
    for i in range(n):
        if (residueId+i) in chain_ends:
            atEnd = True
    return atEnd
def isChainEdge(residueId, chain_starts, chain_ends, n=2):
    # n is how far away from the two ends count as in chain edge.
    return (isChainStart(residueId, chain_starts, n) or isChainEnd(residueId, chain_ends, n))
    # atBegin = False
    # atEnd = False
    # for i in range(n):
    #     if (residueId-i) in chain_starts:
    #         atBegin = True
    # for i in range(n):
    #     if (residueId+i) in chain_ends:
    #         atEnd = True
    # return (atBegin or atEnd)

def inWhichChain(residueId, chain_ends):
    chain_table = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    for i, end_of_chain_resId in enumerate(chain_ends):
        if end_of_chain_resId < residueId:
            pass
        else:
            return chain_table[i]

def read_beta_parameters():
    ### directly copied from Nick Schafer's
    # os.chdir(parameter_directory)
    in_anti_HB = open("anti_HB", 'r').readlines()
    in_anti_NHB = open("anti_NHB", 'r').readlines()
    in_para_HB = open("para_HB", 'r').readlines()
    in_para_one = open("para_one", 'r').readlines()
    in_anti_one = open("anti_one", 'r').readlines()

    p_par = np.zeros((20))
    p_anti = np.zeros((20))
    p_antihb = np.zeros((20,20,2))
    p_antinhb = np.zeros((20,20,2))
    p_parhb = np.zeros((20,20,2))

    for i in range(20):
        p_par[i] = float(in_para_one[i].strip())
        p_anti[i] = float(in_anti_one[i].strip())
        for j in range(20):
            p_antihb[i][j][0] = float(in_anti_HB[i].strip().split()[j])
            p_antinhb[i][j][0] = float(in_anti_NHB[i].strip().split()[j])
            p_parhb[i][j][0] = float(in_para_HB[i].strip().split()[j])

    for i in range(20):
        for j in range(20):
            p_antihb[i][j][1] = float(in_anti_HB[i+21].strip().split()[j])
            p_antinhb[i][j][1] = float(in_anti_NHB[i+21].strip().split()[j])
            p_parhb[i][j][1] = float(in_para_HB[i+21].strip().split()[j])
    return p_par, p_anti, p_antihb, p_antinhb, p_parhb


def get_lambda_by_index(i, j, lambda_i):


    lambda_table = [[1.37, 1.36, 1.17],
                    [3.89, 3.50, 3.52],
                    [0.00, 3.47, 3.62]]
    if abs(j-i) >= 4 and abs(j-i) < 18:
        return lambda_table[lambda_i][0]
    elif abs(j-i) >= 18 and abs(j-i) < 45:
        return lambda_table[lambda_i][1]
    elif abs(j-i) >= 45:
        return lambda_table[lambda_i][2]
    else:
        return 0

def get_alpha_by_index(i, j, alpha_i):
    alpha_table = [[1.30, 1.30, 1.30],
                    [1.32, 1.32, 1.32],
                    [1.22, 1.22, 1.22],
                    [0.00, 0.33, 0.33],
                    [0.00, 1.01, 1.01]]
    if abs(j-i) >= 4 and abs(j-i) < 18:
        return alpha_table[alpha_i][0]
    elif abs(j-i) >= 18 and abs(j-i) < 45:
        return alpha_table[alpha_i][1]
    elif abs(j-i) >= 45:
        return alpha_table[alpha_i][2]
    else:
        return 0

def get_pap_gamma_APH(donor_idx, acceptor_idx, chain_i, chain_j, gamma_APH):
    # if chain_i == chain_j and abs(j-i) < 13 or abs(j-i) > 16:
    # if abs(j-i) < 13 or abs(j-i) > 16:
    # if i-j < 13 or i-j > 16:
    # if (donor_idx - acceptor_idx >= 13 and donor_idx - acceptor_idx <= 16) or chain_i != chain_j:
    if (donor_idx - acceptor_idx >= 13 and donor_idx - acceptor_idx <= 16) and chain_i == chain_j:
        return gamma_APH
    else:
        return 0

def get_pap_gamma_AP(donor_idx, acceptor_idx, chain_i, chain_j, gamma_AP, ssweight):
    if ssweight[donor_idx][1] == 1 and ssweight[acceptor_idx][1] == 1:
        additional_scale = 1.5
    else:
        additional_scale = 1.0
    # if (donor_idx - acceptor_idx >= 17):
    if (donor_idx - acceptor_idx >= 17) or chain_i != chain_j:
        return additional_scale * gamma_AP
    else:
        return 0

def get_pap_gamma_P(donor_idx, acceptor_idx, chain_i, chain_j, gamma_P, ssweight):
    if ssweight[donor_idx][1] == 1 and ssweight[acceptor_idx][1] == 1:
        additional_scale = 1.5
    else:
        additional_scale = 1.0
    if (donor_idx - acceptor_idx >= 9) or chain_i != chain_j:
        return additional_scale * gamma_P
    else:
        return 0

def get_Lambda_2(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a):
    Lambda = get_lambda_by_index(i, j, 1)
    Lambda += -0.5*get_alpha_by_index(i, j, 0)*p_antihb[a[i], a[j]][0]
    Lambda += -0.25*get_alpha_by_index(i, j, 1)*(p_antinhb[a[i+1], a[j-1]][0] + p_antinhb[a[i-1], a[j+1]][0])
    Lambda += -get_alpha_by_index(i, j, 2)*(p_anti[a[i]] + p_anti[a[j]])
    return Lambda

def get_Lambda_3(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a):
    Lambda = get_lambda_by_index(i, j, 2)
    Lambda += -get_alpha_by_index(i, j, 3)*p_parhb[a[i+1], a[j]][0]
    Lambda += -get_alpha_by_index(i, j, 4)*p_par[a[i+1]]
    Lambda += -get_alpha_by_index(i, j, 3)*p_par[a[j]]
    return Lambda


# def beta_term_1(oa, k_beta=4.184):
#     print("beta_1 term ON")
#     nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
#     # print(lambda_1)
#     r_ON = .298
#     sigma_NO = .068
#     r_OH = .206
#     sigma_HO = .076

#     lambda_1 = np.zeros((nres, nres))
#     for i in range(nres):
#         for j in range(nres):
#             lambda_1[i][j] = get_lambda_by_index(i, j, 0)
#     theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
#     beta_string_1 = f"-{k_beta}*lambda_1(res_i,res_j)*theta_ij;theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);"
#     beta_1 = CustomHbondForce(beta_string_1)
#     beta_1.addPerDonorParameter("res_i")
#     beta_1.addPerAcceptorParameter("res_j")
#     beta_1.addTabulatedFunction("lambda_1", Discrete2DFunction(nres, nres, lambda_1.T.flatten()))
#     # print(lambda_1)
#     # print(len(oa.o), nres)
#     for i in range(nres):
#         if oa.o[i]!= -1:
#             beta_1.addAcceptor(oa.o[i], -1, -1, [i])
#         if oa.n[i]!=-1 and oa.h[i]!=-1:
#             beta_1.addDonor(oa.n[i], oa.h[i], -1, [i])
#     beta_1.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#     beta_1.setCutoffDistance(1.0)
#     beta_1.setForceGroup(23)
#     # beta_2.setForceGroup(24)
#     # beta_3.setForceGroup(25)
#     return beta_1

def convert_units(k):
    if isinstance(k, float) or isinstance(k, int):
        k = k   # just for backward comptable
    elif isinstance(k, Quantity):
        k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    else:
        print(f"Unknown input, {k}, {type(k)}")
    return k

def beta_term_1(oa, k=0.5*kilocalories_per_mole, forceGroup=27):
    print("beta_1 term ON")
    k_beta = convert_units(k) * oa.k_awsem
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076

    lambda_1 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            lambda_1[i][j] = get_lambda_by_index(i, j, 0)
    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    mu_1 = 10  # nm^-1
    # mu_2 = 5   # nm^-1
    rcHB = 1.2  # in nm
    # v1i ensures the hydrogen bonding does not occur when five residue segment is shorter than 12 A
    # v1i = f"0.5*(1+tanh({mu_1}*(distance(a2,a3)-{rcHB})))"
    v1i = "1"
    beta_string_1 = f"-{k_beta}*lambda_1(res_i,res_j)*theta_ij*v1i;theta_ij={theta_ij};v1i={v1i};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);"
    beta_1 = CustomHbondForce(beta_string_1)
    beta_1.addPerDonorParameter("res_i")
    beta_1.addPerAcceptorParameter("res_j")
    beta_1.addTabulatedFunction("lambda_1", Discrete2DFunction(nres, nres, lambda_1.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if oa.o[i]!= -1:
            ca_i_minus_2 = oa.ca[0] if i <= 2 else oa.ca[i-2]
            ca_i_plus_2 = oa.ca[-1] if i+2 >= nres else oa.ca[i+2]
            # beta_1.addAcceptor(oa.o[i], ca_i_minus_2, ca_i_plus_2, [i])
            beta_1.addAcceptor(oa.o[i], -1, -1, [i])
        if oa.n[i]!=-1 and oa.h[i]!=-1:
            beta_1.addDonor(oa.n[i], oa.h[i], -1, [i])
    beta_1.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
    beta_1.setCutoffDistance(1.0)
    beta_1.setForceGroup(forceGroup)
    # beta_2.setForceGroup(24)
    # beta_3.setForceGroup(25)
    return beta_1

def beta_term_2(oa, k=0.5*kilocalories_per_mole, forceGroup=27):
    print("beta_2 term ON")
    k_beta = convert_units(k) * oa.k_awsem
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    # for lookup table.
    a = []
    for ii in range(oa.nres):
        a.append(se_map_1_letter[oa.seq[ii]])

    lambda_2 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=1) or \
                    isChainEdge(j, oa.chain_starts, oa.chain_ends, n=1):
                continue
            lambda_2[i][j] = get_Lambda_2(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a)
    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    theta_ji = f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hi-{r_OH})^2/(2*{sigma_HO}^2))"
    beta_string_2 = f"-{k_beta}*lambda_2(res_i,res_j)*theta_ij*theta_ji;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);\
                        theta_ji={theta_ji};r_Oj_Ni=distance(d3,a2);r_Oj_Hi=distance(d3,a3);"
    beta_2 = CustomHbondForce(beta_string_2)
    beta_2.addPerDonorParameter("res_i")
    beta_2.addPerAcceptorParameter("res_j")
    beta_2.addTabulatedFunction("lambda_2", Discrete2DFunction(nres, nres, lambda_2.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if o[i]!= -1 and n[i]!=-1 and h[i]!=-1:
            beta_2.addAcceptor(o[i], n[i], h[i], [i])
            beta_2.addDonor(n[i], h[i], o[i], [i])
    beta_2.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
    beta_2.setCutoffDistance(1.0)
    # beta_1.setForceGroup(23)
    beta_2.setForceGroup(forceGroup)
    # beta_3.setForceGroup(25)

    return beta_2


def beta_term_3(oa, k=0.5*kilocalories_per_mole, forceGroup=27):
    print("beta_3 term ON")
    k_beta = convert_units(k) * oa.k_awsem
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    # for lookup table.
    a = []
    for ii in range(oa.nres):
        a.append(se_map_1_letter[oa.seq[ii]])

    lambda_3 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=1) or \
                    isChainEdge(j, oa.chain_starts, oa.chain_ends, n=1):
                continue
            lambda_3[i][j] = get_Lambda_3(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a)

    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    theta_jip2 = f"exp(-(r_Oj_Nip2-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hip2-{r_OH})^2/(2*{sigma_HO}^2))"

    beta_string_3 = f"-{k_beta}*lambda_3(res_i,res_j)*theta_ij*theta_jip2;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);\
                        theta_jip2={theta_jip2};r_Oj_Nip2=distance(d3,a2);r_Oj_Hip2=distance(d3,a3);"
    beta_3 = CustomHbondForce(beta_string_3)

    beta_3.addPerDonorParameter("res_i")
    beta_3.addPerAcceptorParameter("res_j")
    beta_3.addTabulatedFunction("lambda_3", Discrete2DFunction(nres, nres, lambda_3.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=2):
            continue
        if o[i] != -1 and n[i+2] !=-1 and h[i+2] !=-1:
            beta_3.addAcceptor(o[i], n[i+2], h[i+2], [i])
        if o[i] != -1 and n[i] !=-1 and h[i] !=-1:
            beta_3.addDonor(n[i], h[i], o[i], [i])
    beta_3.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
    beta_3.setCutoffDistance(1.0)
    # beta_1.setForceGroup(23)
    # beta_2.setForceGroup(24)
    beta_3.setForceGroup(forceGroup)

    return beta_3


def pap_term_1(oa, k=0.5*kilocalories_per_mole, dis_i_to_i4=1.2, forceGroup=28, ssweightFileName="ssweight"):
    print("pap_1 term ON")
    k_pap = convert_units(k) * oa.k_awsem
    # dis_i_to_i4 should be in nm, it disfavor hydrogen bond when ca_i and ca_i+4 are 1.2 nm apart away.
    nres, ca = oa.nres, oa.ca
    # r0 = 2.0 # nm
    r0 = 0.8  # nm
    eta_pap = 70  # nm^-1
    gamma_aph = 1.0
    gamma_ap = 0.4
    gamma_p = 0.4

    if not os.path.exists(ssweightFileName):
        print("No ssweight given, assume all zero")
        ssweight = np.zeros((nres, 2))
    else:
        ssweight = np.loadtxt(ssweightFileName)

    gamma_1 = np.zeros((nres, nres))
    gamma_2 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            gamma_1[i][j] = get_pap_gamma_APH(i, j, chain1, chain2, gamma_aph)
            gamma_2[i][j] = get_pap_gamma_AP(i, j, chain1, chain2, gamma_ap, ssweight)

    constraint_i_and_i4 = f"0.5*(1+tanh({eta_pap}*(distance(a1,a2)-{dis_i_to_i4})))"

    pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx)+gamma_2(donor_idx,acceptor_idx))\
                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))\
                        *{constraint_i_and_i4}"

    # pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx))\
    #                     *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
    #                     *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

    # pap_function = f"-{k_pap}*distance(a1,d1)"
    pap = CustomHbondForce(pap_function)
    pap.addPerDonorParameter("donor_idx")
    pap.addPerAcceptorParameter("acceptor_idx")
    pap.addTabulatedFunction("gamma_1", Discrete2DFunction(nres, nres, gamma_1.T.flatten()))
    pap.addTabulatedFunction("gamma_2", Discrete2DFunction(nres, nres, gamma_2.T.flatten()))
    # print(ca)
    # count = 0;
    i = 0

    for i in range(nres):
        if not isChainEnd(i, oa.chain_ends, n=4):
            pap.addAcceptor(ca[i], ca[i+4], -1, [i])

        if not isChainStart(i, oa.chain_starts, n=4):
            if oa.n[i] != -1 and oa.n[i-4] != -1:
                pap.addDonor(oa.n[i], oa.n[i-4], -1, [i])

    pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
    pap.setCutoffDistance(1.0)
    # print(count)
    pap.setForceGroup(forceGroup)
    return pap

def pap_term_2(oa, k=0.5*kilocalories_per_mole, dis_i_to_i4=1.2, forceGroup=28, ssweightFileName="ssweight"):
    print("pap_2 term ON")
    k_pap = convert_units(k) * oa.k_awsem
    nres, ca = oa.nres, oa.ca
    # r0 = 2.0 # nm
    r0 = 0.8  # nm
    eta_pap = 70  # nm^-1
    gamma_aph = 1.0
    gamma_ap = 0.4
    gamma_p = 0.4
    if not os.path.exists(ssweightFileName):
        print("No ssweight given, assume all zero")
        ssweight = np.zeros((nres, 2))
    else:
        ssweight = np.loadtxt(ssweightFileName)

    gamma_3 = np.zeros((nres, nres))
    for i in range(nres):
        for j in range(nres):
            resId1 = i
            chain1 = inWhichChain(resId1, oa.chain_ends)
            resId2 = j
            chain2 = inWhichChain(resId2, oa.chain_ends)
            gamma_3[i][j] = get_pap_gamma_P(i, j, chain1, chain2, gamma_p, ssweight)


    constraint_i_and_i4 = f"0.5*(1+tanh({eta_pap}*(distance(a1,a2)-{dis_i_to_i4})))"
    pap_function = f"-{k_pap}*gamma_3(donor_idx,acceptor_idx)\
                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
                        *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))\
                        *{constraint_i_and_i4}"
    pap = CustomHbondForce(pap_function)
    pap.addPerDonorParameter("donor_idx")
    pap.addPerAcceptorParameter("acceptor_idx")
    pap.addTabulatedFunction("gamma_3", Discrete2DFunction(nres, nres, gamma_3.T.flatten()))
    # print(oa.n)
    # count = 0;
    for i in range(nres):
        if not isChainEnd(i, oa.chain_ends, n=4):
            pap.addAcceptor(ca[i], ca[i+4], -1, [i])
            # pap.addDonor(ca[i], ca[i+4], -1, [i])
            if oa.n[i] != -1 and oa.n[i+4] != -1:
                pap.addDonor(oa.n[i], oa.n[i+4], -1, [i])

    pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
    pap.setCutoffDistance(1.0)
    # print(count)
    pap.setForceGroup(forceGroup)
    return pap

def get_helical_f(oneLetterCode, inMembrane=False):
    if inMembrane:
        table = {"A": 0.79, "R": 0.62, "N": 0.49, "D": 0.44, "C": 0.76, "Q": 0.61, "E": 0.57, "G": 0.57, "H": 0.63, "I": 0.81,
            "L": 0.81, "K": 0.56, "M": 0.80, "F": 0.76, "P": 0.44, "S": 0.6, "T": 0.67, "W": 0.74, "Y": 0.71, "V": 0.79}
    else:
        table = {"A": 0.77, "R": 0.68, "N": 0.07, "D": 0.15, "C": 0.23, "Q": 0.33, "E": 0.27, "G": 0.0, "H": 0.06, "I": 0.23,
            "L": 0.62, "K": 0.65, "M": 0.5, "F": 0.41, "P": 0.4, "S": 0.35, "T": 0.11, "W": 0.45, "Y": 0.17, "V": 0.14}
    return table[oneLetterCode]

def helical_term(oa, k_helical=4.184, inMembrane=False, forceGroup=29):
    # without density dependency.
    # without z dependency for now.
    k_helical *= oa.k_awsem
    sigma_NO = 0.068
    sigma_HO = 0.076
    r_ON = 0.298
    r_OH = 0.206

    theta_ij = f"exp(-(r_Oi_Nip4-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hip4-{r_OH})^2/(2*{sigma_HO}^2))"
    helical = CustomCompoundBondForce(3, f"-{k_helical}*(fa_i+fa_ip4)*{theta_ij};\
                                        r_Oi_Nip4=distance(p1,p2);r_Oi_Hip4=distance(p1,p3);")
    helical.addPerBondParameter("fa_i")
    helical.addPerBondParameter("fa_ip4")
    for i in range(oa.nres):
        # if not isChainEnd(i, oa.chain_ends, n=4) and oa.res_type[i+4] == "IPR":
        #     print(oa.o[i], oa.n[i+4], oa.h[i+4])
        if not isChainEnd(i, oa.chain_ends, n=4) and oa.res_type[i+4] != "IPR":
            fa_i = get_helical_f(oa.seq[i], inMembrane=inMembrane)
            fa_ip4 = get_helical_f(oa.seq[i+4], inMembrane=inMembrane)
            helical.addBond([oa.o[i], oa.n[i+4], oa.h[i+4]], [fa_i, fa_ip4])

    helical.setForceGroup(forceGroup)
    return helical

def z_dependent_helical_term(oa, k_helical=4.184, membrane_center=0*angstrom, z_m=1.5, forceGroup=29):
    # without density dependency.
    k_helical *= oa.k_awsem
    sigma_NO = 0.068
    sigma_HO = 0.076
    r_ON = 0.298
    r_OH = 0.206
    eta_switching = 10
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm

    alpha_membrane = f"0.5*tanh({eta_switching}*((z4-{membrane_center})+{z_m}))+0.5*tanh({eta_switching}*({z_m}-(z4-{membrane_center})))"
    theta_ij = f"exp(-(r_Oi_Nip4-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hip4-{r_OH})^2/(2*{sigma_HO}^2))"
    helical = CustomCompoundBondForce(4, f"-{k_helical}*{theta_ij}*((fa_i+fa_ip4)*(1-alpha_membrane)+(fa_i_membrane+fa_ip4_membrane)*(alpha_membrane));\
                                        alpha_membrane={alpha_membrane};\
                                        r_Oi_Nip4=distance(p1,p2);r_Oi_Hip4=distance(p1,p3);")
    helical.addPerBondParameter("fa_i")
    helical.addPerBondParameter("fa_ip4")
    helical.addPerBondParameter("fa_i_membrane")
    helical.addPerBondParameter("fa_ip4_membrane")
    for i in range(oa.nres):
        # if not isChainEnd(i, oa.chain_ends, n=4) and oa.res_type[i+4] == "IPR":
        #     print(oa.o[i], oa.n[i+4], oa.h[i+4])
        if not isChainEnd(i, oa.chain_ends, n=4) and oa.res_type[i+4] != "IPR":
            fa_i = get_helical_f(oa.seq[i], inMembrane=False)
            fa_ip4 = get_helical_f(oa.seq[i+4], inMembrane=False)
            fa_i_membrane = get_helical_f(oa.seq[i], inMembrane=True)
            fa_ip4_membrane = get_helical_f(oa.seq[i+4], inMembrane=True)
            helical.addBond([oa.o[i], oa.n[i+4], oa.h[i+4], oa.ca[i]], [fa_i, fa_ip4, fa_i_membrane, fa_ip4_membrane])

    helical.setForceGroup(forceGroup)
    return helical


# def pap_term_1(oa, k_pap=4.184, dis_i_to_i4=-1):
#     print("pap_1 term ON")
#     nres, ca = oa.nres, oa.ca
#     # r0 = 2.0 # nm
#     r0 = 0.8 # nm
#     eta_pap = 70 # nm^-1
#     gamma_aph = 1.0
#     gamma_ap = 0.4
#     gamma_p = 0.4

#     gamma_1 = np.zeros((nres, nres))
#     gamma_2 = np.zeros((nres, nres))
#     for i in range(nres):
#         for j in range(nres):
#             resId1 = i
#             chain1 = inWhichChain(resId1, oa.chain_ends)
#             resId2 = j
#             chain2 = inWhichChain(resId2, oa.chain_ends)
#             gamma_1[i][j] = get_pap_gamma_APH(i, j, chain1, chain2, gamma_aph)
#             gamma_2[i][j] = get_pap_gamma_AP(i, j, chain1, chain2, gamma_ap)

#     pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx)+gamma_2(donor_idx,acceptor_idx))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

#     # pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx))\
#     #                     *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#     #                     *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

#     # pap_function = f"-{k_pap}*distance(a1,d1)"
#     pap = CustomHbondForce(pap_function)
#     pap.addPerDonorParameter("donor_idx")
#     pap.addPerAcceptorParameter("acceptor_idx")
#     pap.addTabulatedFunction("gamma_1", Discrete2DFunction(nres, nres, gamma_1.T.flatten()))
#     pap.addTabulatedFunction("gamma_2", Discrete2DFunction(nres, nres, gamma_2.T.flatten()))
#     # print(ca)
#     # count = 0;
#     i = 0
#     # pap.addAcceptor(ca[0], ca[4], -1, [0])
#     # pap.addAcceptor(ca[20], ca[8], -1, [4])
#     # pap.addDonor(ca[20], ca[0], -1, [4])
#     for i in range(nres):
#         if not isChainEnd(i, oa.chain_ends, n=4):
#             pap.addAcceptor(ca[i], ca[i+4], -1, [i])

#         if i > 13 and not isChainStart(i, oa.chain_starts, n=4):
#             pap.addDonor(oa.n[i], oa.n[i-4], -1, [i])

#     pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#     pap.setCutoffDistance(1.0)
#     # print(count)
#     pap.setForceGroup(26)
#     return pap


# def pap_term_1(oa, k_pap=4.184):
#     print("pap_1 term ON")
#     nres, ca = oa.nres, oa.ca
#     # r0 = 2.0 # nm
#     r0 = 0.8 # nm
#     eta_pap = 70 # nm^-1
#     gamma_aph = 1.0
#     gamma_ap = 0.4
#     gamma_p = 0.4

#     gamma_1 = np.zeros((nres, nres))
#     gamma_2 = np.zeros((nres, nres))
#     for i in range(nres):
#         for j in range(nres):
#             resId1 = i
#             chain1 = inWhichChain(resId1, oa.chain_ends)
#             resId2 = j
#             chain2 = inWhichChain(resId2, oa.chain_ends)
#             gamma_1[i][j] = get_pap_gamma_APH(i, j, chain1, chain2, gamma_aph)
#             gamma_2[i][j] = get_pap_gamma_AP(i, j, chain1, chain2, gamma_ap)

#     pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx)+gamma_2(donor_idx,acceptor_idx))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

#     pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

#     pap_function = f"-{k_pap}*distance(a1,d1)"
#     pap = CustomHbondForce(pap_function)
#     pap.addPerDonorParameter("donor_idx")
#     pap.addPerAcceptorParameter("acceptor_idx")
#     pap.addTabulatedFunction("gamma_1", Discrete2DFunction(nres, nres, gamma_1.T.flatten()))
#     pap.addTabulatedFunction("gamma_2", Discrete2DFunction(nres, nres, gamma_2.T.flatten()))
#     # print(ca)
#     # count = 0;
#     i = 0
#     # pap.addAcceptor(ca[0], ca[4], -1, [0])
#     # pap.addAcceptor(ca[20], ca[8], -1, [4])
#     # pap.addDonor(ca[20], ca[0], -1, [4])
#     for i in range(nres):
#         if not isChainEnd(i, oa.chain_ends, n=4):
#             pap.addAcceptor(ca[i], ca[i+4], -1, [i])

#         if i > 13 and not isChainStart(i, oa.chain_starts, n=4):
#             pap.addDonor(ca[i], ca[i-4], -1, [i])

#     pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#     pap.setCutoffDistance(1.0)
#     # print(count)
#     pap.setForceGroup(26)
#     return pap


def beta_term_1_old(oa, k_beta=4.184, debug=False, forceGroup=23):

    print("beta_1 term ON")
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # print(lambda_1)
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2

    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    # theta_ji =   f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hi-{r_OH})^2/(2*{sigma_HO}^2))"
    # theta_jip2 = "exp(-(r_Oj_Nip2-r_ON)^2/(2*sigma_NO^2)-(r_Oj_Hip2-r_OH)^2/(2*sigma_HO^2))"
    nu_i = f"0.5*(1+tanh({eta_beta_1}*(r_CAim2_CAip2-{r_HB_c})))"
    nu_j = f"0.5*(1+tanh({eta_beta_2}*(r_CAjm2_CAjp2-{r_HB_c})))"

    # Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4     5     6     7
    beta_string_1 = f"-k_beta*lambda_1*theta_ij*nu_i*nu_j;theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                    nu_i={nu_i};nu_j={nu_j};r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)"
    # # below used for debug, set, vi vj = 0
    if debug:
        beta_string_1 = f"-k_beta*lambda_1*theta_ij*nu_i*nu_j;theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                        nu_i=1+0*{nu_i};nu_j=1+0*{nu_j};r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)"

    # beta_string_1 = f"-k_beta*lambda_1"
    # beta_string_1 = f"-k_beta"

    beta_1 = CustomCompoundBondForce(7, beta_string_1)
    # beta_2 = CustomCompoundBondForce(10, beta_string_2)
    # beta_3 = CustomCompoundBondForce(10, beta_string_3)
    # add parameters to force
    beta_1.addGlobalParameter("k_beta", k_beta)
    beta_1.addPerBondParameter("lambda_1")
    # beta_2.addTabulatedFunction("lambda_2", Discrete2DFunction(nres, nres, lambda_2))
    # beta_3.addTabulatedFunction("lambda_3", Discrete2DFunction(nres, nres, lambda_3))

    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=2) or \
                isChainEdge(j, oa.chain_starts, oa.chain_ends, n=2):
                continue
            if not res_type[j] == "IPR":
                beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [get_lambda_by_index(i, j, 0)])
            #if not res_type[i] == "IPR" and not res_type[j] == "IPR":
            #    beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
            #if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
            #    beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])

    # beta_1.setForceGroup(23)
    #beta_2.setForceGroup(24)
    #beta_3.setForceGroup(25)
    beta_1.setForceGroup(forceGroup)
    return beta_1

def beta_term_2_old(oa, k_beta=4.184, debug=False, forceGroup=24):
    print("beta_2 term ON");
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # add beta potential
    # setup parameters
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    theta_ij =   f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    theta_ji =   f"exp(-(r_Oj_Ni-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hi-{r_OH})^2/(2*{sigma_HO}^2))"
    nu_i = f"0.5*(1+tanh({eta_beta_1}*(r_CAim2_CAip2-{r_HB_c})))"
    nu_j = f"0.5*(1+tanh({eta_beta_2}*(r_CAjm2_CAjp2-{r_HB_c})))"

    # Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4     5     6     7
    #beta_string_1 = "-k_beta*lambda_1(index_i,index_j)*theta_ij*nu_i*nu_j;theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
    #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)" % (theta_ij, nu_i, nu_j)

    # Oi Nj Hj Oj Ni Hi CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4  5  6  7     8     9     10
    beta_string_2 = f"-k_beta*lambda_2*theta_ij*theta_ji*nu_i*nu_j;\
                    theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                    theta_ji={theta_ji};r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
                    nu_i={nu_i};nu_j={nu_j};r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)"
    # # below used for debug, set, vi vj = 0
    if debug:
        beta_string_2 = f"-k_beta*lambda_2*theta_ij*theta_ji*nu_i*nu_j;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                        theta_ji={theta_ji};r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
                        nu_i=1+0*{nu_i};nu_j=1+0*{nu_j};r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)"

    # Oi Nj Hj Oj Ni+2 Hi+2 CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4  5    6    7     8     9     10
    #beta_string_3 = "-k_beta*lambda_3(index_i,index_j)*theta_ij*theta_jip2*nu_i*nu_j;\
    #                theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
    #                theta_ji=%s;r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
    #                theta_jip2=%s;r_Oj_Nip2=distance(p4,p5);r_Oj_Hip2=distance(p4,p6);\
    #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)" % (theta_ij, theta_ji, theta_jip2, nu_i, nu_j)

    #beta_1 = CustomCompoundBondForce(7, beta_string_1)
    beta_2 = CustomCompoundBondForce(10, beta_string_2)
    #beta_3 = CustomCompoundBondForce(10, beta_string_3)
    # add parameters to force
    beta_2.addGlobalParameter("k_beta", k_beta)
    beta_2.addPerBondParameter("lambda_2")

    # for lookup table.
    a = []
    for ii in range(oa.nres):
        a.append(se_map_1_letter[oa.seq[ii]])

    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=2) or \
                isChainEdge(j, oa.chain_starts, oa.chain_ends, n=2):
                continue
            #if not res_type[j] == "IPR":
            #    beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
            if not res_type[i] == "IPR" and not res_type[j] == "IPR":
                beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [get_Lambda_2(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a)])
            #if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
            #    beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])


    #beta_1.setForceGroup(23)
    beta_2.setForceGroup(forceGroup)
    #beta_3.setForceGroup(25)
    return beta_2

def beta_term_3_old(oa, k_beta=4.184, debug=False, forceGroup=25):
    print("beta_3 term ON")
    nres, n, h, ca, o, res_type = oa.nres, oa.n, oa.h, oa.ca, oa.o, oa.res_type
    # add beta potential
    # setup parameters
    r_ON = .298
    sigma_NO = .068
    r_OH = .206
    sigma_HO = .076
    eta_beta_1 = 10.0
    eta_beta_2 = 5.0
    # r_HB_c = 0.4
    r_HB_c = 1.2
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()

    theta_ij = f"exp(-(r_Oi_Nj-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oi_Hj-{r_OH})^2/(2*{sigma_HO}^2))"
    theta_jip2 = f"exp(-(r_Oj_Nip2-{r_ON})^2/(2*{sigma_NO}^2)-(r_Oj_Hip2-{r_OH})^2/(2*{sigma_HO}^2))"
    nu_i = f"0.5*(1+tanh({eta_beta_1}*(r_CAim2_CAip2-{r_HB_c})))"
    nu_j = f"0.5*(1+tanh({eta_beta_2}*(r_CAjm2_CAjp2-{r_HB_c})))"

    # Oi Nj Hj CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4     5     6     7
    #beta_string_1 = "-k_beta*lambda_1(index_i,index_j)*theta_ij*nu_i*nu_j;theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
    #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p4,p5);r_CAjm2_CAjp2=distance(p6,p7)" % (theta_ij, nu_i, nu_j)

    # Oi Nj Hj Oj Ni Hi CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4  5  6  7     8     9     10
    #beta_string_2 = "-k_beta*lambda_2(index_i,index_j)*theta_ij*theta_ji*nu_i*nu_j;\
    #                theta_ij=%s;r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
    #                theta_ji=%s;r_Oj_Ni=distance(p4,p5);r_Oj_Hi=distance(p4,p6);\
    #                nu_i=%s;nu_j=%s;r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)" % (theta_ij, theta_ji, nu_i, nu_j)

    # Oi Nj Hj Oj Ni+2 Hi+2 CAi-2 CAi+2 CAj-2 CAj+2
    # 1  2  3  4  5    6    7     8     9     10
    beta_string_3 = f"-k_beta*lambda_3*theta_ij*theta_jip2*nu_i*nu_j;\
                    theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                    theta_jip2={theta_jip2};r_Oj_Nip2=distance(p4,p5);r_Oj_Hip2=distance(p4,p6);\
                    nu_i={nu_i};nu_j={nu_j};r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)"
    # # below used for debug, set, vi vj = 0
    if debug:
        beta_string_3 = f"-k_beta*lambda_3*theta_ij*theta_jip2*nu_i*nu_j;\
                        theta_ij={theta_ij};r_Oi_Nj=distance(p1,p2);r_Oi_Hj=distance(p1,p3);\
                        theta_jip2={theta_jip2};r_Oj_Nip2=distance(p4,p5);r_Oj_Hip2=distance(p4,p6);\
                        nu_i=1+0*{nu_i};nu_j=1+0*{nu_j};r_CAim2_CAip2=distance(p7,p8);r_CAjm2_CAjp2=distance(p9,p10)"

    beta_3 = CustomCompoundBondForce(10, beta_string_3)
    # add parameters to force
    beta_3.addGlobalParameter("k_beta", k_beta)
    beta_3.addPerBondParameter("lambda_3")

    # for lookup table.
    a = []
    for ii in range(oa.nres):
        a.append(se_map_1_letter[oa.seq[ii]])

    for i in range(nres):
        for j in range(nres):
            if isChainEdge(i, oa.chain_starts, oa.chain_ends, n=2) or \
                isChainEdge(j, oa.chain_starts, oa.chain_ends, n=2):
                continue
            #if not res_type[j] == "IPR":
            #    beta_1.addBond([o[i], n[j], h[j], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
            #if not res_type[i] == "IPR" and not res_type[j] == "IPR":
            #    beta_2.addBond([o[i], n[j], h[j], o[j], n[i], h[i], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [i, j])
            if not res_type[i+2] == "IPR" and not res_type[j] == "IPR":
                beta_3.addBond([o[i], n[j], h[j], o[j], n[i+2], h[i+2], ca[i-2], ca[i+2], ca[j-2], ca[j+2]], [get_Lambda_3(i, j, p_par, p_anti, p_antihb, p_antinhb, p_parhb, a)])


    #beta_1.setForceGroup(23)
    #beta_2.setForceGroup(24)
    beta_3.setForceGroup(forceGroup)
    return beta_3

def pap_term_old(oa, k_pap=4.184, forceGroup=26):
    print("pap term ON")
    nres, ca = oa.nres, oa.ca
    # r0 = 2.0 # nm
    r0 = 0.8 # nm
    eta_pap = 70 # nm^-1
    gamma_aph = 1.0
    gamma_ap = 0.4
    gamma_p = 0.4
    pap_function = f"-k_pap*gamma*0.5*(1+tanh({eta_pap}*({r0}-distance(p1,p2))))*0.5*(1+tanh({eta_pap}*({r0}-distance(p3,p4))))"
    pap = CustomCompoundBondForce(4, pap_function)
    pap.addGlobalParameter("k_pap", k_pap)
    pap.addPerBondParameter("gamma")
    #count = 0;
    for i in range(nres):
        for j in range(nres):
            # anti-parallel hairpin for i from 1 to N-13 and j from i+13 to min(i+16,N)
            # CAi CAj CAi+4 CAj-4
            # 1   2   3     4
            if i <= nres-13 and j >= i+13 and j <= min(i+16,nres):
                pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_aph])
                #count = count + 1
                #print([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_aph])
            # anti-parallel for i from 1 to N-17 and j from i+17 to N
            # CAi CAj CAi+4 CAj-4
            # 1   2   3     4
            if i <= nres-17 and j >= i+17 and j <= nres:
                pap.addBond([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_ap])
                #count = count + 1;
                #print([ca[i], ca[j], ca[i+4], ca[j-4]], [gamma_ap])
            # parallel for i from 1 to N-13 and j from i+9 to N-4
            # CAi CAj CAi+4 CAj+4
            # 1   2   3     4
            if i <= nres-13 and j >= i+9 and j < nres-4:
                #print([i, j, i+4, j+4])
                #print([i, j, i+4, j+4, ca[i], ca[j], ca[i+4], ca[j+4]], [gamma_p])
                pap.addBond([ca[i], ca[j], ca[i+4], ca[j+4]], [gamma_p])
                #count = count + 1;

    # print(count)
    pap.setForceGroup(forceGroup)
    return pap

# def pap_term_1(oa, k_pap=4.184):
#     print("pap_1 term ON")
#     nres, ca = oa.nres, oa.ca
#     # r0 = 2.0 # nm
#     r0 = 0.8 # nm
#     eta_pap = 70 # nm^-1
#     gamma_aph = 1.0
#     gamma_ap = 0.4
#     gamma_p = 0.4

#     gamma_1 = np.zeros((nres, nres))
#     gamma_2 = np.zeros((nres, nres))
#     for i in range(nres):
#         for j in range(nres):
#             resId1 = i
#             chain1 = inWhichChain(resId1, oa.chain_ends)
#             resId2 = j
#             chain2 = inWhichChain(resId2, oa.chain_ends)
#             gamma_1[i][j] = get_pap_gamma_APH(i, j, chain1, chain2, gamma_aph)
#             gamma_2[i][j] = get_pap_gamma_AP(i, j, chain1, chain2, gamma_ap)

#     pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx)+gamma_2(donor_idx,acceptor_idx))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

#     pap_function = f"-{k_pap}*(gamma_1(donor_idx,acceptor_idx))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a1,d1))))\
#                         *0.5*(1+tanh({eta_pap}*({r0}-distance(a2,d2))))"

#     pap_function = f"-{k_pap}*distance(a1,d1)"
#     pap = CustomHbondForce(pap_function)
#     pap.addPerDonorParameter("donor_idx")
#     pap.addPerAcceptorParameter("acceptor_idx")
#     pap.addTabulatedFunction("gamma_1", Discrete2DFunction(nres, nres, gamma_1.T.flatten()))
#     pap.addTabulatedFunction("gamma_2", Discrete2DFunction(nres, nres, gamma_2.T.flatten()))
#     # print(ca)
#     # count = 0;
#     i = 0
#     # pap.addAcceptor(ca[0], ca[4], -1, [0])
#     # pap.addAcceptor(ca[20], ca[8], -1, [4])
#     # pap.addDonor(ca[20], ca[0], -1, [4])
#     for i in range(nres):
#         if not isChainEnd(i, oa.chain_ends, n=4):
#             pap.addAcceptor(ca[i], ca[i+4], -1, [i])

#         if i > 13 and not isChainStart(i, oa.chain_starts, n=4):
#             pap.addDonor(ca[i], ca[i-4], -1, [i])

#     pap.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
#     pap.setCutoffDistance(1.0)
#     # print(count)
#     pap.setForceGroup(26)
#     return pap


'''
# old way of getting lambda
def lambda_coefficient(oa, i, j, lambda_index):
    p_par, p_anti, p_antihb, p_antinhb, p_parhb = read_beta_parameters()
    parameter_i = []
    # print(i,j,lambda_index)
    for ii in range(oa.nres):
        # print(oa.seq[i])
        parameter_i.append(se_map_1_letter[oa.seq[ii]])
    # print(p_antihb[parameter_i[i], parameter_i[j]][0],p_antinhb[parameter_i[i+1],parameter_i[j-1]][0],p_anti[parameter_i[i]], p_anti[parameter_i[j]])
    lambda_2_extra_terms = -0.5*oa.alpha_coefficient(parameter_i[i],parameter_i[j],1)*p_antihb[parameter_i[i], parameter_i[j]][0]-0.25*oa.alpha_coefficient(parameter_i[i], parameter_i[j], 2)*(p_antinhb[parameter_i[i+1],parameter_i[j-1]][0] + p_antinhb[parameter_i[i-1],parameter_i[j+1]][0])-oa.alpha_coefficient(parameter_i[i], parameter_i[j], 3)*(p_anti[parameter_i[i]]+p_anti[parameter_i[j]])
    lambda_3_extra_terms = -oa.alpha_coefficient(parameter_i[i],parameter_i[j], 4)*p_parhb[parameter_i[i+1],parameter_i[j]][0]-oa.alpha_coefficient(parameter_i[i],parameter_i[j],5)*p_par[parameter_i[i+1]]+oa.alpha_coefficient(parameter_i[i],parameter_i[j],4)*p_par[parameter_i[j]]
    if abs(j-i) >= 4 and abs(j-i) < 18:
        if lambda_index == 1:
            return 1.37
        elif lambda_index == 2:
            return 3.89+lambda_2_extra_terms
        elif lambda_index == 3:
            return 0.0+lambda_3_extra_terms
    elif abs(j-i) >= 18 and abs(j-i) < 45:
        if lambda_index == 1:
            return 1.36
        elif lambda_index == 2:
            return 3.50+lambda_2_extra_terms
        elif lambda_index == 3:
            return 3.47+lambda_3_extra_terms
    elif abs(j-i) >= 45:
        if lambda_index == 1:
            return 1.17
        elif lambda_index == 2:
            return 3.52+lambda_2_extra_terms
        elif lambda_index == 3:
            return 3.62+lambda_3_extra_terms
    elif abs(j-i) < 4:
        return 0.0

def alpha_coefficient(oa, i,j, alpha_index):
    if abs(j-i) >= 4 and abs(j-i) < 18:
        if alpha_index == 1:
            return 1.3
        if alpha_index == 2:
            return 1.32
        if alpha_index == 3:
            return 1.22
        if alpha_index == 4:
            return 0.0
        if alpha_index == 5:
            return 0.0
    elif abs(j-i) >= 18 and abs(j-i) < 45:
        if alpha_index == 1:
            return 1.3
        if alpha_index == 2:
            return 1.32
        if alpha_index == 3:
            return 1.22
        if alpha_index == 4:
            return 0.33
        if alpha_index == 5:
            return 1.01
    elif abs(j-i) >= 45:
        if alpha_index == 1:
            return 1.3
        if alpha_index == 2:
            return 1.32
        if alpha_index == 3:
            return 1.22
        if alpha_index == 4:
            return 0.33
        if alpha_index == 5:
            return 1.01
    elif abs(j-i) <4:
        return 0.0
'''
