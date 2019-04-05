from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def membrane_term(oa, k_membrane=4.184, k_m=2, z_m=1.5):
    # k_m in units of nm^-1, z_m in units of nm.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k_membrane *= oa.k_awsem
    membrane = CustomExternalForce(f"{k_membrane}*\
            (0.5*tanh({k_m}*(z+{z_m}))+0.5*tanh({k_m}*({z_m}-z)))*hydrophobicityScale")
    membrane.addPerParticleParameter("hydrophobicityScale")
    zim = np.loadtxt("zim")
    cb_fixed = [x if x > 0 else y for x,y in zip(oa.cb,oa.ca)]
    for i in cb_fixed:
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
        membrane.addParticle(i, [zim[oa.resi[i]]])
    membrane.setForceGroup(20)
    return membrane

