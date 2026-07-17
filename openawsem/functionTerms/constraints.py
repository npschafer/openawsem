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

def measure_distance(oa, res1, res2, forceGroup=4): #Assign to forceGroup 4 as measurement placeholder; Rg measurement is RESERVED forceGroup 3.
    # print(len(oa.ca))
    constraint = CustomBondForce(f"(r)")
    # res1, res2 is 0 index. res1 = 0 means the first residue.
    constraint.addBond(*[oa.ca[res1], oa.ca[res2]])         # you could also do constraint.addBond(oa.ca[res1], oa.ca[res2])
    constraint.setForceGroup(forceGroup)
    return constraint


def group_constraint_by_distance(oa, d0=0*angstrom, group1=None, group2=None, forceGroup=3, k=1*kilocalorie_per_mole):
    # example: group1 = [oa.ca[0], oa.ca[1]], group2 = [oa.ca[2], oa.ca[3]]
    # CustomCentroidBondForce only work with CUDA not OpenCL.
    #
    # note added 11 Jun 2025: CustomCentroidBondForce worked for me on OpenCL on my workstation, ws1808
    #
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


def group_constraint_by_distance_protein(oa, d0=0*angstrom, group1=None, group2=None, forceGroup=3, k=1*kilocalorie_per_mole):
    # CustomCentroidBondForce only work with CUDA not OpenCL.
    # only CA, CB, O has mass. so the group have to include those. Default assignment should be away from forceGroup 3.
    if group1 is None or group2 is None:
        raise ValueError("Both group1 and group2 must be provided as lists of residue indices.")
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    d0 = d0.value_in_unit(nanometer)   # convert to nm
    constraint = CustomCentroidBondForce(2, f"0.5*{k_constraint}*(distance(g1,g2)-{d0})^2")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    residues1 = []
    residues2 = []
    for r in group1:
        for a in [oa.ca[r], oa.cb[r], oa.o[r]]:
        #for a in [oa.ca[r]]:
            if a != -1:    # Catch glycine CB cases
                residues1.append(a)
                #print(f"Added residue index {r} and atom index {a}")
    for r in group2:
        for a in [oa.ca[r], oa.cb[r], oa.o[r]]:
        #for a in [oa.ca[r]]:
            if a != -1:    # Catch glycine CB cases
                residues2.append(a)
                #print(f"Added residue index {r} and atom index {a}")
    #print(f"Group 1 initially has {len(group1)} atoms and {len(residues1)} in.")
    #print(f"Group 1 initially has {len(group2)} atoms and {len(residues2)} in.")
    constraint.addGroup(residues1)    # group use particle index.
    constraint.addGroup(residues2)
    constraint.addBond([0, 1])
    constraint.setForceGroup(forceGroup)
    return constraint

def measure_distance_group(oa, group1=None, group2=None, forceGroup=4): #Assign to forceGroup 4 as measurement placeholder; Rg measurement is RESERVED forceGroup 3.
    if group1 is None or group2 is None:
        raise ValueError("Both group1 and group2 must be provided as lists of residue indices.")
    residues1 = []
    residues2 = []
    for r in group1:
        for a in [oa.ca[r], oa.cb[r], oa.o[r]]:
        #for a in [oa.ca[r]]:
            if a != -1:    # Catch glycine CB cases
                residues1.append(a)
                #print(f"Added residue index {r} and atom index {a}")
    for r in group2:
        for a in [oa.ca[r], oa.cb[r], oa.o[r]]:
        #for a in [oa.ca[r]]:
            if a != -1:    # Catch glycine CB cases
                residues2.append(a)
                #print(f"Added residue index {r} and atom index {a}")
    constraint = CustomCentroidBondForce(2, f"distance(g1,g2)")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    constraint.addGroup(residues1)    # group use particle index.
    constraint.addGroup(residues2)
    constraint.addBond([0, 1])
    constraint.setForceGroup(forceGroup)
    return constraint

def group_constraint_by_position(oa, k=1*kilocalorie_per_mole, x0=10*nanometer, y0=10*nanometer, z0=10*nanometer, appliedToResidues=None, forceGroup=3):
    # x0, y0, z0 is in unit of nm.
    x0 = x0.value_in_unit(nanometer)
    y0 = y0.value_in_unit(nanometer)
    z0 = z0.value_in_unit(nanometer)
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
        if appliedToResidues is None:
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

def group_constraint_by_position_centroid(oa,k,restraint_coords,restraint_group,forceGroup=30):
    """
    k: kcal/mol, force constant for the restraint, where the potential is (1/2)*k*x^2
    restraint_coords: nm, a list of coordinates [x,y,z] to which the group is restrained
    restraint_group: a list of 0-indexed residue indices of the group to be restrained 
    forceGroup: i think default value of 30 is not used by any other openawsem term, but am not sure
    """
    # initialize Force
    restraint_force = CustomCentroidBondForce(1,f"0.5*{k}*sqrt(((x1-{restraint_coords[0]})^2)+((y1-{restraint_coords[1]})^2)+((z1-{restraint_coords[2]})^2))")
    # get particle indices
    particles = []
    for counter,ca_i in enumerate(oa.ca):
        if ca_i != -1 and counter in restraint_group:
            particles.append(ca_i)
    for counter,cb_i in enumerate(oa.cb):
        if cb_i != -1 and counter in restraint_group:
            particles.append(cb_i)
    for counter,o_i in enumerate(oa.o):
        if o_i != -1 and counter in restraint_group:
            particles.append(o_i)
    # add particles to Force
    restraint_force.addGroup(particles)
    #
    restraint_force.setForceGroup(forceGroup)
    return restraint_force


def measure_from_position(oa, x0=10*angstrom, y0=10*angstrom, z0=10*angstrom, appliedToResidues=None, forceGroup=3):
    # x0, y0, z0 is in unit of nm.
    x0 = x0.value_in_unit(nanometer)
    y0 = y0.value_in_unit(nanometer)
    z0 = z0.value_in_unit(nanometer)
    # appliedToResidues can be a list of residue index. for example appliedToResidues=[0, 1], to tether the first two residues.
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
        if appliedToResidues == None:
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
    harmonic = CustomCVForce(f"(sum_x/{total_mass}-{x0})^2+(sum_y/{total_mass}-{y0})^2+(sum_z/{total_mass}-{z0})^2")
    harmonic.addCollectiveVariable("sum_x", sum_of_x_coord)
    harmonic.addCollectiveVariable("sum_y", sum_of_y_coord)
    harmonic.addCollectiveVariable("sum_z", sum_of_z_coord)
    harmonic.setForceGroup(forceGroup)
    return harmonic

# For more advanced OpenMM users, you can try using atom/coarse grained particle index.
def group_index_constraint_by_distance(oa, d0=0*angstrom, group1=None, group2=None, forceGroup=3, k=1*kilocalorie_per_mole):
    if group1 is None or group2 is None:
        raise ValueError("Both group1 and group2 must be provided as lists of particle indices.")
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_constraint = k * oa.k_awsem
    d0 = d0.value_in_unit(nanometer)   # convert to nm
    constraint = CustomCentroidBondForce(2, f"0.5*{k_constraint}*(distance(g1,g2)-{d0})^2")
    # example group set up group1=[oa.ca[7], oa.cb[7]] use the ca and cb of residue 8.
    #print(f"Group 1 initially has {len(group1)} atoms and {len(residues1)} in.")
    #print(f"Group 1 initially has {len(group2)} atoms and {len(residues2)} in.")
    constraint.addGroup(group1)    # group use particle index.
    constraint.addGroup(group2)
    constraint.addBond([0, 1])
    constraint.setForceGroup(forceGroup)
    return constraint

def measure_distance_group_index(oa, group1=None, group2=None, forceGroup=4): #Assign to forceGroup 4 as measurement placeholder; Rg measurement is RESERVED forceGroup 3.
    #oa argument is not used at all here but is a placeholder to simplify forces setup file and avoid errors namely double assignment.
    #The first argument in every openawsem function is the oa function which takes in the protein itself.
    if group1 is None or group2 is None:
        raise ValueError("Both group1 and group2 must be provided as lists of particle indices.")
    constraint = CustomCentroidBondForce(2, f"distance(g1,g2)")
    constraint.addGroup(group1)    # group use particle index.
    constraint.addGroup(group2)
    constraint.addBond([0, 1])
    constraint.setForceGroup(forceGroup)
    return constraint

def group_index_constraint_by_position(oa, k=1*kilocalorie_per_mole, x0=10*angstrom, y0=10*angstrom, z0=10*angstrom, appliedToResidues=None, forceGroup=3):
    # appliedToResidues really takes in Particle Incidies and NOT Residue Indicies; variable name left as is for now.
    # x0, y0, z0 is in unit of nm.
    x0 = x0.value_in_unit(nanometer)
    y0 = y0.value_in_unit(nanometer)
    z0 = z0.value_in_unit(nanometer)
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
    #print(f"mass can be retrieved as ", oa.system.getParticleMass(oa.ca[0]))
    total_mass = 0.0
    if appliedToResidues == None:
        for i in range(oa.natoms):
            #mass = oa.system.getParticleMass(i).value_in_unit(dalton)
            mass = 1   #mass = 1 is a temporary placeholder
            sum_of_x_coord.addParticle(i, [mass])
            sum_of_y_coord.addParticle(i, [mass])
            sum_of_z_coord.addParticle(i, [mass])
            total_mass += mass
    else:
        for i in appliedToResidues:
            #mass = oa.system.getParticleMass(i).value_in_unit(dalton)
            mass = 1
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

def measure_from_position_index(oa, x0=10*angstrom, y0=10*angstrom, z0=10*angstrom, appliedToResidues=None, forceGroup=4):
    # appliedToResidues really takes in Particle Incidies and NOT Residue Indicies; variable name left as is for now.
    # x0, y0, z0 is in unit of nm.
    x0 = x0.value_in_unit(nanometer)
    y0 = y0.value_in_unit(nanometer)
    z0 = z0.value_in_unit(nanometer)
    # appliedToResidues can be a list of residue index. for example appliedToResidues=[0, 1], to tether the first two residues.
    sum_of_x_coord = CustomExternalForce(f"x*mass")
    sum_of_y_coord = CustomExternalForce(f"y*mass")
    sum_of_z_coord = CustomExternalForce(f"z*mass")

    sum_of_x_coord.addPerParticleParameter("mass")
    sum_of_y_coord.addPerParticleParameter("mass")
    sum_of_z_coord.addPerParticleParameter("mass")

    # print("index for CAs", oa.ca)
    #print(f"mass can be retrieved as ", oa.system.getParticleMass(oa.ca[0]))
    total_mass = 0.0
    if appliedToResidues == None:
        for i in range(oa.natoms):
            #mass = oa.system.getParticleMass(i).value_in_unit(dalton)
            mass = 1
            sum_of_x_coord.addParticle(i, [mass])
            sum_of_y_coord.addParticle(i, [mass])
            sum_of_z_coord.addParticle(i, [mass])
            total_mass += mass
    else:
        for i in appliedToResidues:
            #mass = oa.system.getParticleMass(i).value_in_unit(dalton)
            mass = 1
            sum_of_x_coord.addParticle(i, [mass])
            sum_of_y_coord.addParticle(i, [mass])
            sum_of_z_coord.addParticle(i, [mass])
            total_mass += mass
        # if oa.resi[i] == appliedToResidue:
        #     pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    print(f"total_mass = {total_mass}")
    harmonic = CustomCVForce(f"(sum_x/{total_mass}-{x0})^2+(sum_y/{total_mass}-{y0})^2+(sum_z/{total_mass}-{z0})^2")
    harmonic.addCollectiveVariable("sum_x", sum_of_x_coord)
    harmonic.addCollectiveVariable("sum_y", sum_of_y_coord)
    harmonic.addCollectiveVariable("sum_z", sum_of_z_coord)
    harmonic.setForceGroup(forceGroup)
    return harmonic

# This function is to allow for openawsem-style force setup of openmm orientational constraints class. Input includes particle indicies, not oa protein residue indicies.
# For now, will constrain the orientation to initial positions of the particles selected.
def orientational_constraints(oa, k = 100*kilocalorie_per_mole, particles = None, forceGroup=5):
    k = k.value_in_unit(kilojoule_per_mole)   # must be converted to kJ/mol

    if particles is None:
        orient_force = OrientationRestraintForce(k, oa.pdb.positions)
    else:
        orient_force = OrientationRestraintForce(k, oa.pdb.positions, particles)
    
    orient_force.setForceGroup(forceGroup)

    return orient_force
