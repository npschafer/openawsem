#!/usr/bin/env python3
import math

# line = "{type:4}  {serialNumber:5} {AtomName:4} {ReidueName:3} {Chain:1}{ResId:4} {x:8.3f}{y:8.3f}{z:8.3f}{1:6.2f}{0:6.2f}    {element} \n"


def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def read_lammps(lammps_file="dump.lammpstrj", center=False, ca=True):
    nFrame = 0
    with open(lammps_file, "r") as lfile:
        for line in lfile:
            l = line.strip()
            if l[:5]=="ITEM:":
                item = l[6:]
            else:
                if item == "TIMESTEP":
                    step = int(l)
                    try:
                        atoms_all_frames.append(atoms)
                    except NameError:
                        atoms_all_frames = []
                    atoms = []
                    box = []
                    A = []
                    nFrame = nFrame + 1
                elif item == "NUMBER OF ATOMS":
                        n_atoms = int(l)
                        xyz_count = 0
                elif item[:10] == "BOX BOUNDS":
                    # I center x, y not z
                    if center:
                        if xyz_count <= 1:
                            xyz_count += 1
                            box.append(l)
                            l = l.split()
                            # A.append([float(l[0]), float(l[1])])
                            l_left = (float(l[0]) - float(l[1]))/2.0
                            l_right = (float(l[1]) - float(l[0]))/2.0
                            A.append([l_left, l_right])
                            # print l_right - l_left
                        else:
                            xyz_count = 0
                            box.append(l)
                            l = l.split()
                            A.append([float(l[0]), float(l[1])])
                            # l_left = (float(l[0]) - float(l[1]))/2.0
                            # l_right = (float(l[1]) - float(l[0]))/2.0
                            # A.append([l_left, l_right])
                            # print l_right - l_left
                    else:
                        box.append(l)
                        l = l.split()
                        A.append([float(l[0]), float(l[1])])
                elif item[:5] == "ATOMS":
                    l = l.split()
                    i_atom = int(l[0])
                    x = float(l[2])
                    y = float(l[3])
                    z = float(l[4])
                    x = (A[0][1] - A[0][0])*x + A[0][0]
                    y = (A[1][1] - A[1][0])*y + A[1][0]
                    z = (A[2][1] - A[2][0])*z + A[2][0]
                    # C alpha distance
                    if ca:
                        if i_atom % 3 == 1:
                            atom = [x, y, z]
                            atoms.append(atom)
                    else:   # C beta or H in the case of GLY
                        if i_atom % 3 == 0:
                            atom = [x, y, z]
                            atoms.append(atom)
        atoms_all_frames.append(atoms)
    return atoms_all_frames
