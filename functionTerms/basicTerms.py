from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def con_term(oa, k_con=50208, bond_lengths=[.3816, .240, .276, .153], forceGroup=20):
    # add con forces
    # 50208 = 60 * 2 * 4.184 * 100. kJ/nm^2, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_con *= oa.k_awsem
    con = HarmonicBondForce()
    for i in range(oa.nres):
        con.addBond(oa.ca[i], oa.o[i], bond_lengths[1], k_con)
        if not oa.res_type[i] == "IGL":  # OpenAWSEM system doesn't have virtual HB atom for glycine, so con_term doesn't include C_alpha and HB interaction in LAMMPS
            con.addBond(oa.ca[i], oa.cb[i], bond_lengths[3], k_con)
        if i not in oa.chain_ends:
            con.addBond(oa.ca[i], oa.ca[i+1], bond_lengths[0], k_con)
            con.addBond(oa.o[i], oa.ca[i+1], bond_lengths[2], k_con)
    con.setForceGroup(forceGroup)   # start with 11, so that first 10 leave for user defined force.
    return con


def chain_term(oa, k_chain=50208, bond_lengths=[0.2459108, 0.2519591, 0.2466597], forceGroup=20):
    # add chain forces
    # 50208 = 60 * 2 * 4.184 * 100. kJ/nm^2, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_chain *= oa.k_awsem
    chain = HarmonicBondForce()
    for i in range(oa.nres):
        if i not in oa.chain_starts and not oa.res_type[i] == "IGL":
            chain.addBond(oa.n[i], oa.cb[i], bond_lengths[0], k_chain)
        if i not in oa.chain_ends and not oa.res_type[i] == "IGL":
            chain.addBond(oa.c[i], oa.cb[i], bond_lengths[1], k_chain)
        if i not in oa.chain_starts and i not in oa.chain_ends:
            chain.addBond(oa.n[i], oa.c[i], bond_lengths[2], k_chain)
    chain.setForceGroup(forceGroup)
    return chain

def chi_term(oa, k_chi=251.04, chi0=-0.71, forceGroup=20):
    # add chi forces
    # The sign of the equilibrium value is opposite and magnitude differs slightly
    # 251.04 = 60 * 4.184 kJ, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_chi *= oa.k_awsem
    chi = CustomCompoundBondForce(4, "k_chi*(chi*norm-chi0)^2;"
                                        "chi=crossproduct_x*r_cacb_x+crossproduct_y*r_cacb_y+crossproduct_z*r_cacb_z;"
                                        "crossproduct_x=(u_y*v_z-u_z*v_y);"
                                        "crossproduct_y=(u_z*v_x-u_x*v_z);"
                                        "crossproduct_z=(u_x*v_y-u_y*v_x);"
                                        "norm=1/((u_x*u_x+u_y*u_y+u_z*u_z)*(v_x*v_x+v_y*v_y+v_z*v_z)*(r_cacb_x*r_cacb_x+r_cacb_y*r_cacb_y+r_cacb_z*r_cacb_z))^0.5;"
                                        "r_cacb_x=x1-x4;"
                                        "r_cacb_y=y1-y4;"
                                        "r_cacb_z=z1-z4;"
                                        "u_x=x1-x2; u_y=y1-y2; u_z=z1-z2;"
                                        "v_x=x3-x1; v_y=y3-y1; v_z=z3-z1;")
    chi.addGlobalParameter("k_chi", k_chi)
    chi.addGlobalParameter("chi0", chi0)
    for i in range(oa.nres):
        if i not in oa.chain_starts and i not in oa.chain_ends and not oa.res_type[i] == "IGL":
            chi.addBond([oa.ca[i], oa.c[i], oa.n[i], oa.cb[i]])
    chi.setForceGroup(forceGroup)
    return chi

def excl_term(oa, k_excl=8368, r_excl=0.35, periodic=False, forceGroup=20):
    # add excluded volume
    # Still need to add element specific parameters
    # 8368 = 20 * 4.184 * 100 kJ/nm^2, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    # Openawsem doesn't have the distance range (r_excl) change from 0.35 to 0.45 when the sequence separtation more than 5
    k_excl *= oa.k_awsem
    excl = CustomNonbondedForce(f"{k_excl}*step({r_excl}-r)*(r-{r_excl})^2")
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

def rama_term(oa, k_rama=8.368, num_rama_wells=3, w=[1.3149, 1.32016, 1.0264], sigma=[15.398, 49.0521, 49.0954], omega_phi=[0.15, 0.25, 0.65], phi_i=[-1.74, -1.265, 1.041], omega_psi=[0.65, 0.45, 0.25], psi_i=[2.138, -0.318, 0.78], forceGroup=21):
    # add Rama potential
    # 8.368 = 2 * 4.184 kJ/mol, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_rama *= oa.k_awsem
    rama_function = ''.join(["w%d*exp(-sigma%d*(omega_phi%d*phi_term%d^2+omega_psi%d*psi_term%d^2))+" \
                            % (i, i, i, i, i, i) for i in range(num_rama_wells)])[:-1]
    rama_function = '-k_rama*(' + rama_function + ");"
    rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                            psi_term{i}=cos(psi_{i}-psi0{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                            for i in range(num_rama_wells)])
    rama_string = rama_function+rama_parameters
    rama = CustomCompoundBondForce(5, rama_string)
    for i in range(num_rama_wells):
        rama.addGlobalParameter(f"k_rama", k_rama)
        rama.addGlobalParameter(f"w{i}", w[i])
        rama.addGlobalParameter(f"sigma{i}", sigma[i])
        rama.addGlobalParameter(f"omega_phi{i}", omega_phi[i])
        rama.addGlobalParameter(f"omega_psi{i}", omega_psi[i])
        rama.addGlobalParameter(f"phi0{i}", phi_i[i])
        rama.addGlobalParameter(f"psi0{i}", psi_i[i])
    for i in range(oa.nres):
        if i not in oa.chain_starts and i not in oa.chain_ends and not oa.res_type[i] == "IGL" and not oa.res_type[i] == "IPR":
            rama.addBond([oa.c[i-1], oa.n[i], oa.ca[i], oa.c[i], oa.n[i+1]])
    rama.setForceGroup(forceGroup)
    return rama

def rama_proline_term(oa, k_rama_proline=8.368, num_rama_proline_wells=2, w=[2.17, 2.15], sigma=[105.52, 109.09], omega_phi=[1.0, 1.0], phi_i=[-1.153, -0.95], omega_psi=[0.15, 0.15], psi_i=[2.4, -0.218], forceGroup=21):
    # add Rama potential for prolines
    # 8.368 = 2 * 4.184 kJ/mol, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_rama_proline *= oa.k_awsem
    rama_function = ''.join(["w_P%d*exp(-sigma_P%d*(omega_phi_P%d*phi_term%d^2+omega_psi_P%d*psi_term%d^2))+" \
                            % (i, i, i, i, i, i) for i in range(num_rama_proline_wells)])[:-1]
    rama_function = '-k_rama_proline*(' + rama_function + ");"
    rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0_P{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                            psi_term{i}=cos(psi_{i}-psi0_P{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                            for i in range(num_rama_proline_wells)])
    rama_string = rama_function+rama_parameters
    rama = CustomCompoundBondForce(5, rama_string)
    for i in range(num_rama_proline_wells):
        rama.addGlobalParameter(f"k_rama_proline", k_rama_proline)
        rama.addGlobalParameter(f"w_P{i}", w[i])
        rama.addGlobalParameter(f"sigma_P{i}", sigma[i])
        rama.addGlobalParameter(f"omega_phi_P{i}", omega_phi[i])
        rama.addGlobalParameter(f"omega_psi_P{i}", omega_psi[i])
        rama.addGlobalParameter(f"phi0_P{i}", phi_i[i])
        rama.addGlobalParameter(f"psi0_P{i}", psi_i[i])
    for i in range(oa.nres):
        if i not in oa.chain_starts and i not in oa.chain_ends and oa.res_type[i] == "IPR":
            rama.addBond([oa.c[i-1], oa.n[i], oa.ca[i], oa.c[i], oa.n[i+1]])
    rama.setForceGroup(forceGroup)
    return rama

def rama_ssweight_term(oa, k_rama_ssweight=8.368, num_rama_wells=2, w=[2.0, 2.0],
                    sigma=[419.0, 15.398], omega_phi=[1.0, 1.0], phi_i=[-0.995, -2.25],
                    omega_psi=[1.0, 1.0], psi_i=[-0.82, 2.16], location_pre="./", forceGroup=21):
    # add RamaSS potential
    # 8.368 = 2 * 4.184 kJ/mol, converted from default value in LAMMPS AWSEM
    # multiply interaction strength by overall scaling
    k_rama_ssweight *= oa.k_awsem
    rama_function = ''.join(["wSS%d*ssweight(%d,resId)*exp(-sigmaSS%d*(omega_phiSS%d*phi_term%d^2+omega_psiSS%d*psi_term%d^2))+" \
                            % (i, i, i, i, i, i, i) for i in range(num_rama_wells)])[:-1]
    rama_function = f'-{k_rama_ssweight}*(' + rama_function + ");"
    rama_parameters = ''.join([f"phi_term{i}=cos(phi_{i}-phi0SS{i})-1; phi_{i}=dihedral(p1, p2, p3, p4);\
                            psi_term{i}=cos(psi_{i}-psi0SS{i})-1; psi_{i}=dihedral(p2, p3, p4, p5);"\
                            for i in range(num_rama_wells)])
    rama_string = rama_function+rama_parameters
    ramaSS = CustomCompoundBondForce(5, rama_string)
    ramaSS.addPerBondParameter("resId")
    for i in range(num_rama_wells):
        ramaSS.addGlobalParameter(f"wSS{i}", w[i])
        ramaSS.addGlobalParameter(f"sigmaSS{i}", sigma[i])
        ramaSS.addGlobalParameter(f"omega_phiSS{i}", omega_phi[i])
        ramaSS.addGlobalParameter(f"omega_psiSS{i}", omega_psi[i])
        ramaSS.addGlobalParameter(f"phi0SS{i}", phi_i[i])
        ramaSS.addGlobalParameter(f"psi0SS{i}", psi_i[i])
    for i in range(oa.nres):
        if i not in oa.chain_starts and i not in oa.chain_ends and not oa.res_type[i] == "IGL" and not oa.res_type == "IPR":
            ramaSS.addBond([oa.c[i-1], oa.n[i], oa.ca[i], oa.c[i], oa.n[i+1]], [i])
    ssweight = np.loadtxt(location_pre+"ssweight")
    ramaSS.addTabulatedFunction("ssweight", Discrete2DFunction(2, oa.nres, ssweight.flatten()))
    ramaSS.setForceGroup(forceGroup)
    return ramaSS
