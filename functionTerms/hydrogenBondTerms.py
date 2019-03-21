from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

se_map_1_letter = {'A': 0,  'P': 1,  'K': 2,  'N': 3,  'R': 4,
                   'F': 5,  'D': 6,  'Q': 7,  'E': 8,  'G': 9,
                   'I': 10, 'H': 11, 'L': 12, 'C': 13, 'M': 14,
                   'S': 15, 'T': 16, 'Y': 17, 'V': 18, 'W': 19}

def isChainEdge(residueId, chain_starts, chain_ends, n=2):
    # n is how far away from the two ends count as in chain edge.
    atBegin = False
    atEnd = False
    for i in range(n):
        if (residueId-i) in chain_starts:
            atBegin = True
    for i in range(n):
        if (residueId+i) in chain_ends:
            atEnd = True
    return (atBegin or atEnd)

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


def beta_term_1(oa, k_beta=4.184):
    print("beta_1 term ON")
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
    beta_string_1 = f"-{k_beta}*lambda_1(res_i,res_j)*theta_ij;theta_ij={theta_ij};r_Oi_Nj=distance(a1,d1);r_Oi_Hj=distance(a1,d2);"
    beta_1 = CustomHbondForce(beta_string_1)
    beta_1.addPerDonorParameter("res_i")
    beta_1.addPerAcceptorParameter("res_j")
    beta_1.addTabulatedFunction("lambda_1", Discrete2DFunction(nres, nres, lambda_1.T.flatten()))
    # print(lambda_1)
    # print(len(oa.o), nres)
    for i in range(nres):
        if oa.o[i]!= -1:
            beta_1.addAcceptor(oa.o[i], -1, -1, [i])
        if oa.n[i]!=-1 and oa.h[i]!=-1:
            beta_1.addDonor(oa.n[i], oa.h[i], -1, [i])
    beta_1.setNonbondedMethod(CustomHbondForce.CutoffNonPeriodic)
    beta_1.setCutoffDistance(1.0)
    beta_1.setForceGroup(23)
    # beta_2.setForceGroup(24)
    # beta_3.setForceGroup(25)
    return beta_1

def beta_term_2(oa, k_beta=4.184):
    print("beta_2 term ON")
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
    beta_2.setForceGroup(24)
    # beta_3.setForceGroup(25)

    return beta_2


def beta_term_3(oa, k_beta=4.184):
    print("beta_3 term ON")
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
    beta_3.setForceGroup(25)

    return beta_3


def beta_term_1_old(oa, k_beta=4.184, debug=False):

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

    beta_1.setForceGroup(23)
    #beta_2.setForceGroup(24)
    #beta_3.setForceGroup(25)
    return beta_1

def beta_term_2_old(oa, k_beta=4.184, debug=False):
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
    beta_2.setForceGroup(24)
    #beta_3.setForceGroup(25)
    return beta_2

def beta_term_3_old(oa, k_beta=4.184, debug=False):
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
    beta_3.setForceGroup(25)
    return beta_3

def pap_term(oa, k_pap=4.184):
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

    #print(count)
    pap.setForceGroup(26)
    return pap

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
