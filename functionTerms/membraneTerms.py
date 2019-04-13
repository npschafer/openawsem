from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def membrane_term(oa, k_membrane=4.184, k_m=2, z_m=1.5, membrane_center=0):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k_membrane *= oa.k_awsem
    membrane = CustomExternalForce(f"{k_membrane}*\
            (0.5*tanh({k_m}*((z-{membrane_center})+{z_m}))+0.5*tanh({k_m}*({z_m}-(z-{membrane_center}))))*hydrophobicityScale")
    membrane.addPerParticleParameter("hydrophobicityScale")
    zim = np.loadtxt("zim")
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    for i in cb_fixed:
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
        membrane.addParticle(i, [zim[oa.resi[i]]])
    membrane.setForceGroup(20)
    return membrane

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
