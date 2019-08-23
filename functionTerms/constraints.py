from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np


def group_constraint_by_distance(oa, group1=[0], group2=[10], forceGroup=3, k=1*kilocalorie_per_mole):
    # CustomCentroidBondForce only work with CUDA not OpenCL.
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    constraint = CustomCentroidBondForce(2, f"0.5*{k_constraint}*distance(g1,g2)^2")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    constraint.addGroup(group1)    # group use particle index.
    constraint.addGroup(group2)
    constraint.addBond([0, 1])
    constraint.setForceGroup(forceGroup)
    return constraint
