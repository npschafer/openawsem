try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
import numpy as np

def constraint_by_distance(oa, res1, res2,  d0=0*angstrom, forceGroup=3, k=1*kilocalorie_per_mole):
    # print(len(oa.ca))
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    d0 = d0.value_in_unit(nanometer)   # convert to nm
    constraint = CustomBondForce(f"0.5*{k_constraint}*(r-{d0})^2")
    # res1, res2 is 0 index. res1 = 0 means the first residue.
    constraint.addBond(*[oa.ca[res1], oa.ca[res2]])         # you could also do constraint.addBond(oa.ca[res1], oa.ca[res2])
    constraint.setForceGroup(forceGroup)
    return constraint

def group_constraint_by_distance(oa, d0=0*angstrom, group1=[oa.ca[0], oa.ca[1]], group2=[oa.ca[2], oa.ca[3]], forceGroup=3, k=1*kilocalorie_per_mole):
    # CustomCentroidBondForce only work with CUDA not OpenCL.
    # only CA, CB, O has mass. so the group have to include those.
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    d0 = d0.value_in_unit(nanometer)   # convert to nm
    constraint = CustomCentroidBondForce(2, f"0.5*{k_constraint}*(distance(g1,g2)-{d0})^2")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    constraint.addGroup(group1)    # group use particle index.
    constraint.addGroup(group2)
    constraint.addBond([0, 1])
    constraint.setForceGroup(forceGroup)
    return constraint

from simtk.unit import dalton
def group_constraint_by_position(oa, k=1*kilocalorie_per_mole, x0=10, y0=10, z0=10, appliedToResidues=-1, forceGroup=3):
    # x0, y0, z0 is in unit of nm.
    # appliedToResidues can be a list of residue index. for example appliedToResidues=[0, 1], to tether the first two residues.
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    sum_of_x_coord = CustomExternalForce(f"x*mass")
    sum_of_y_coord = CustomExternalForce(f"y*mass")
    sum_of_z_coord = CustomExternalForce(f"z*mass")

    sum_of_x_coord.addPerParticleParameter("mass")
    sum_of_y_coord.addPerParticleParameter("mass")
    sum_of_z_coord.addPerParticleParameter("mass")

    # print("index for CAs", oa.ca)
    print(f"mass can be retrieved as ", oa.system.getParticleMass(oa.ca[0]))
    total_mass = 0.0
    for i in range(oa.natoms):
        if appliedToResidues == -1:
            mass = oa.system.getParticleMass(i).value_in_unit(dalton)
            sum_of_x_coord.addParticle(i, [mass])
            sum_of_y_coord.addParticle(i, [mass])
            sum_of_z_coord.addParticle(i, [mass])
            total_mass += mass
        elif oa.resi[i] in appliedToResidues:
            mass = oa.system.getParticleMass(i).value_in_unit(dalton)
            sum_of_x_coord.addParticle(i, [mass])
            sum_of_y_coord.addParticle(i, [mass])
            sum_of_z_coord.addParticle(i, [mass])
            total_mass += mass
        # if oa.resi[i] == appliedToResidue:
        #     pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    print(f"total_mass = {total_mass}")
    harmonic = CustomCVForce(f"{k_constraint}*((sum_x/{total_mass}-{x0})^2+(sum_y/{total_mass}-{y0})^2+(sum_z/{total_mass}-{z0})^2)")
    harmonic.addCollectiveVariable("sum_x", sum_of_x_coord)
    harmonic.addCollectiveVariable("sum_y", sum_of_y_coord)
    harmonic.addCollectiveVariable("sum_z", sum_of_z_coord)
    harmonic.setForceGroup(forceGroup)
    return harmonic

