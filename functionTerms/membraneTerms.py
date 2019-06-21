from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def membrane_term(oa, k=1*kilocalorie_per_mole, k_m=20, z_m=1.5, membrane_center=0*angstrom):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling

    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_membrane = k * oa.k_awsem

    membrane = CustomExternalForce(f"k_membrane*\
            (0.5*tanh({k_m}*((z-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z-{membrane_center}))))*hydrophobicityScale-0.5")
    membrane.addPerParticleParameter("hydrophobicityScale")
    membrane.addGlobalParameter("k_membrane", k_membrane)
    zim = np.loadtxt("zim")
    # cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    ca = oa.ca
    for i in ca:
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
        membrane.addParticle(i, [zim[oa.resi[i]]])
    membrane.setForceGroup(20)
    return membrane



def membrane_preassigned_term(oa, k_membrane=4.184, k_m=20, z_m=1.5, membrane_center=0):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k_membrane *= oa.k_awsem
    membrane = CustomExternalForce(f"{k_membrane}*\
            (0.5*tanh({k_m}*((z-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z-{membrane_center}))))*zim-0.5")
    membrane.addPerParticleParameter("zim")
    # zim = np.loadtxt("zim")
    zimPosition = np.loadtxt("zimPosition")
    zim = [-1 if z == 2 else 1 for z in zimPosition]
    # print(zim)
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    for i in cb_fixed:
        membrane.addParticle(i, [zim[oa.resi[i]]])
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    membrane.setForceGroup(20)
    return membrane



def single_helix_orientation_bias_term(oa, k=1*kilocalorie_per_mole, membrane_center=0*angstrom, z_m=1.5, k_m=20, atomGroup=-1):
    membrane_center = membrane_center.value_in_unit(nanometer)   # convert to nm
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_single_helix_orientation_bias = oa.k_awsem * k
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup
    n = len(group)
    theta_z1 = f"(0.5*tanh({k_m}*((z1-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z1-{membrane_center}))))"
    theta_z2 = f"(0.5*tanh({k_m}*((z2-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z2-{membrane_center}))))"
    normalization = n * n
    v_orientation = CustomCompoundBondForce(2, f"k_single_helix_orientation_bias/{normalization}*((x1-x2)^2+(y1-y2)^2)*{theta_z1}*{theta_z2}")
    # rcm_square = CustomCompoundBondForce(2, "1/normalization*(x1*x2)")
    v_orientation.addGlobalParameter("k_single_helix_orientation_bias", k_single_helix_orientation_bias)
    # rg_square = CustomBondForce("1/normalization*(sqrt(x^2+y^2)-rcm))^2")
    # rg = CustomBondForce("1")
    # v_orientation.addGlobalParameter("normalization", n*n)
    for i in group:
        for j in group:
            if j <= i:
                continue
            v_orientation.addBond([ca[i], ca[j]], [])

    v_orientation.setForceGroup(27)
    return v_orientation


def pulling_term(oa, k_pulling=4.184, forceDirect="x", appliedToResidue=0):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k_pulling *= oa.k_awsem
    pulling = CustomExternalForce(f"(-{k_pulling})*({forceDirect})")
    for i in range(oa.natoms):
        if oa.resi[i] == appliedToResidue:
            pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    pulling.setForceGroup(29)
    return pulling

def rg_term(oa, convertToAngstrom=True):
    rg_square = CustomBondForce("1/normalization*r^2")
    # rg = CustomBondForce("1")
    rg_square.addGlobalParameter("normalization", oa.nres*oa.nres)
    for i in range(oa.nres):
        for j in range(i+1, oa.nres):
            rg_square.addBond(oa.ca[i], oa.ca[j], [])
    if convertToAngstrom:
        unit = 10
    else:
        unit = 1
    rg = CustomCVForce(f"{unit}*rg_square^0.5")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(2)
    return rg

def rg_bias_term(oa, k_rg=4.184, rg0=0, atomGroup=-1):
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup
    n = len(group)
    rg_square = CustomBondForce("1/normalization*r^2")
    # rg = CustomBondForce("1")
    rg_square.addGlobalParameter("normalization", n*n)
    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond(ca[i], ca[j], [])
    rg = CustomCVForce(f"{k_rg}*(rg_square^0.5-{rg0})^2")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(27)
    return rg