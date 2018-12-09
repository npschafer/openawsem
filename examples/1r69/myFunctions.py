#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import subprocess
import glob
import re
from helperFunctions.myFunctions_helper import *
import numpy as np
import pandas as pd
import fileinput
from itertools import product
from Bio.PDB.PDBParser import PDBParser

from Bio.PDB import PDBList
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

# compute cross Q for every pdb pair in one folder
# parser = argparse.ArgumentParser(description="Compute cross q")
# parser.add_argument("-m", "--mode",
#                     type=int, default=1)

# args = parser.parse_args()

def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()

def read_hydrophobicity_scale(seq, isNew=False):
    seq_dataFrame = pd.DataFrame({"oneLetterCode":list(seq)})
    HFscales = pd.read_table("~/opt/small_script/Whole_residue_HFscales.txt")
    if not isNew:
        # Octanol Scale
        # new and old difference is at HIS.
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS+" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    else:
        code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
                "ARG+" : "R", "LYS+" : "K", "MET" : "M", "CYS" : "C",
                "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
                "TRP" : "W", "ASP-" : "D", "GLU-" : "E", "ASN" : "N",
                "GLN" : "Q", "PHE" : "F", "HIS0" : "H", "VAL" : "V",
                "M3L" : "K", "MSE" : "M", "CAS" : "C"}
    HFscales_with_oneLetterCode = HFscales.assign(oneLetterCode=HFscales.AA.str.upper().map(code)).dropna()
    data = seq_dataFrame.merge(HFscales_with_oneLetterCode, on="oneLetterCode", how="left")
    return data

def create_zim(seqFile, isNew=False):
    a = seqFile
    seq = getFromTerminal("cat " + a).rstrip()
    data = read_hydrophobicity_scale(seq, isNew=isNew)
    z = data["DGwoct"].values
    np.savetxt("zim", z, fmt="%.2f")


def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())


def duplicate_pdb(From, To, offset_x=0, offset_y=0, offset_z=0, new_chain="B"):
    with open(To, "w") as out:
        with open(From, "r") as f:
            for line in f:
                tmp = list(line)
                atom = line[0:4]
                atomSerialNumber = line[6:11]
                atomName = line[12:16]
                atomResidueName = line[17:20]
                chain = line[21]
                residueNumber = line[22:26]
                # change chain A to B
                # new_chain = "B"
                tmp[21] = new_chain
                if atom == "ATOM":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # add 40 to the x
                    new_x = x + offset_x
                    new_y = y + offset_y
                    new_z = z + offset_z


                    tmp[30:38] = "{:8.3f}".format(new_x)
                    tmp[38:46] = "{:8.3f}".format(new_y)
                    tmp[46:54] = "{:8.3f}".format(new_z)

                a = "".join(tmp)
                out.write(a)

def compute_native_contacts(coords, MAX_OFFSET=4, DISTANCE_CUTOFF=9.5):
    native_coords = np.array(coords)
    a= native_coords[:,np.newaxis]
    dis = np.sqrt(np.sum((a - native_coords)**2, axis=2))

    n = len(dis)
    remove_band = np.eye(n)
    for i in range(1, MAX_OFFSET):
        remove_band += np.eye(n, k=i)
        remove_band += np.eye(n, k=-i)
    dis[remove_band==1] = np.max(dis)
    native_contacts = dis < DISTANCE_CUTOFF
    return native_contacts.astype("int")

def compute_contacts(coords, native_contacts, DISTANCE_CUTOFF=9.5):
    native_coords = np.array(coords)
    a= native_coords[:,np.newaxis]
    dis = np.sqrt(np.sum((a - native_coords)**2, axis=2))
    constacts = dis < DISTANCE_CUTOFF
    constacts = constacts*native_contacts  # remove non native contacts
    return np.sum(constacts, axis=1).astype("float")

def compute_localQ_init(MAX_OFFSET=4, DISTANCE_CUTOFF=9.5):
    from pathlib import Path
    home = str(Path.home())
    struct_id = '2xov'
    filename = os.path.join(home, "opt/pulling/2xov.pdb")
    p = PDBParser(PERMISSIVE=1)
    s = p.get_structure(struct_id, filename)
    chains = s[0].get_list()

    # import pdb file
    native_coords = []
    for chain in chains:
        dis = []
        all_res = []
        for res in chain:
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            if (res.get_resname()=='GLY'):
                native_coords.append(res['CA'].get_coord())
            elif (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
                native_coords.append(res['CB'].get_coord())
            else:
                print('ERROR: irregular residue at %s!' % res)
                exit()
    native_contacts_table = compute_native_contacts(native_coords, MAX_OFFSET, DISTANCE_CUTOFF)

    return native_contacts_table

def compute_localQ(native_contacts_table, pre=".", ii=-1, MAX_OFFSET=4, DISTANCE_CUTOFF=9.5):
    native_contacts = np.sum(native_contacts_table, axis=1).astype("float")
    dump = read_lammps(os.path.join(pre, f"dump.lammpstrj.{ii}"), ca=False)
    localQ_list = []
    for atom in dump:
        contacts = compute_contacts(np.array(atom), native_contacts_table, DISTANCE_CUTOFF=DISTANCE_CUTOFF)
        c = np.divide(contacts, native_contacts, out=np.zeros_like(contacts), where=native_contacts!=0)
        localQ_list.append(c)
    data = pd.DataFrame(localQ_list)
    data.columns = ["Res" + str(i+1) for i in data.columns]
    data.to_csv(os.path.join(pre, f"localQ.{ii}.csv"), index=False)

def readPMF_basic(pre):
    # perturbation_table = {0:"original", 1:"p_mem",
    #                       2:"m_mem", 3:"p_lipid",
    #                       4:"m_lipid", 5:"p_go",
    #                       6:"m_go", 7:"p_rg", 8:"m_rg"}
    perturbation_table = {0:"original", 1:"m_go",
                          2:"p_go", 3:"m_lipid",
                          4:"p_lipid", 5:"m_mem",
                          6:"p_mem", 7:"m_rg", 8:"p_rg"}
    pmf_list = {
        "perturbation":list(perturbation_table.keys())
    }
    pmf_list_data = expand_grid(pmf_list)
    all_pmf_list = []
    for index, row in pmf_list_data.iterrows():
        perturbation = row["perturbation"]
        if perturbation == 0:
            location = pre + f"/pmf-*.dat"
            pmf_list = glob.glob(location)
            change = "none"
            upOrDown = "none"
        else:
            location = pre + f"/perturbation-{perturbation}-pmf-*.dat"
            pmf_list = glob.glob(location)
            change = perturbation_table[perturbation].split("_")[-1]
            upOrDown = perturbation_table[perturbation].split("_")[0]
        # print(location)
        name_list = ["f", "df", "e", "s"]
        names = ["bin", "x"] + name_list
        for location in pmf_list:
            # print(location)
            temp = re.findall(r'pmf-(\d+)', location)
            if len(temp) != 1:
                raise ValueError('Not expected to see more than one or none')
            else:
                temp = temp[0]
            data = pd.read_table(location, skiprows=2, sep='\s+', names=names).assign(upOrDown=upOrDown, change=change, temp=temp, perturbation=perturbation_table[perturbation])
            all_pmf_list.append(data)

    return pd.concat(all_pmf_list).dropna().reset_index()

def make_metadata_3(k=1000.0, temps_list=["450"], i=-1, biasLow=None, biasHigh=None):
    print("make metadata")
    cwd = os.getcwd()
    files = glob.glob(f"../data_{i}/*")
    kconstant = k
    with open("metadatafile", "w") as out:
        for oneFile in sorted(files):
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            bias = tmp.split("_")[3]
            if biasLow:
                if float(bias) < biasLow:
                    continue
            if biasHigh:
                if float(bias) > biasHigh:
                    continue
            # print(tmp)
            # if int(float(dis)) > 150:
            #     continue
            if t in temps_list:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, bias)
                out.write(target)


def readPMF(pre, is2d=False, force_list=["0.0", "0.1", "0.2"]):
    # perturbation_table = {0:"original", 1:"p_mem",
    #                       2:"m_mem", 3:"p_lipid",
    #                       4:"m_lipid", 5:"p_go",
    #                       6:"m_go", 7:"p_rg", 8:"m_rg"}
    perturbation_table = {0:"original", 1:"m_go",
                          2:"p_go", 3:"m_lipid",
                          4:"p_lipid", 5:"m_mem",
                          6:"p_mem", 7:"m_rg", 8:"p_rg"}
    pmf_list = {
        "perturbation":list(perturbation_table.keys()),
        "force":force_list
    }
    pmf_list_data = expand_grid(pmf_list)
    all_pmf_list = []
    for index, row in pmf_list_data.iterrows():
        force = row["force"]
        perturbation = row["perturbation"]
        if perturbation == 0:
            location = pre + f"/force_{force}/pmf-*.dat"
            pmf_list = glob.glob(location)
            change = "none"
            upOrDown = "none"
        else:
            location = pre + f"/force_{force}/perturbation-{perturbation}-pmf-*.dat"
            pmf_list = glob.glob(location)
            change = perturbation_table[perturbation].split("_")[-1]
            upOrDown = perturbation_table[perturbation].split("_")[0]
        # print(pmf_list)
        name_list = ["f", "df", "e", "s"]
        if is2d:
            names = ["x", "y"] + name_list
        else:
            names = ["bin", "x"] + name_list
        for location in pmf_list:
            # print(location)
            temp = re.findall(r'pmf-(\d+)', location)
            if len(temp) != 1:
                raise ValueError('Not expected to see more than one or none')
            else:
                temp = temp[0]
            data = pd.read_table(location, skiprows=2, sep='\s+', names=names).assign(upOrDown=upOrDown, change=change, force=force, temp=temp, perturbation=perturbation_table[perturbation])
            all_pmf_list.append(data)

    return pd.concat(all_pmf_list).dropna().reset_index()

def readPMF_2(pre, is2d=0, force_list=["0.0", "0.1", "0.2"]):
    if is2d:
        print("reading 2d pmfs")
    else:
        print("reading 1d dis, qw and z")
    if is2d == 1:
        mode_list = ["2d_qw_dis", "2d_z_dis", "2d_z_qw"]
    elif is2d == 2:
        mode_list = ["quick"]
    else:
        mode_list = ["1d_dis", "1d_qw", "1d_z"]
    all_data_list =[]
    for mode in mode_list:
        tmp = readPMF(mode, is2d, force_list).assign(mode=mode)
        all_data_list.append(tmp)
    return pd.concat(all_data_list).dropna().reset_index()

def shrinkage(n=552, shrink_size=6, max_frame=2000, fileName="dump.lammpstrj"):
    print("Shrinkage: size: {}, max_frame: {}".format(shrink_size, max_frame))
    bashCommand = "wc " + fileName
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    line_number = int(output.decode("utf-8").split()[0])
    print(line_number)
    print(line_number/552)
    # number of atom = 543
    n = 552
    count = 0
    with open("small.lammpstrj", "w") as out:
        with open(fileName, "r") as f:
            for i, line in enumerate(f):
                if (i // n) % shrink_size == 0:
                    if count >= max_frame*n:
                        break
                    count += 1
                    out.write(line)

def compute_theta_for_each_helix(output="angles.csv", dumpName="../dump.lammpstrj.0"):
    print("This is for 2xov only")
    helices_list = [(94,114), (147,168), (171, 192), (200, 217), (226, 241), (250, 269)]
    atoms_all_frames = read_lammps(dumpName)
    # print(atoms[0])
    # print(len(atoms), len(atoms[0]))
    # helices_angles_all_frames = []
    with open(output, "w") as out:
        out.write("Frame, Helix, Angle\n")
        for ii, frame in enumerate(atoms_all_frames):
            # helices_angles = []
            for count, (i, j) in enumerate(helices_list):
                # print(i, j)
                i = i-91
                j = j-91
                # end - start
                a = np.array(frame[j]) - np.array(frame[i])
                b = np.array([0, 0, 1])
                angle = a[2]/length(a)  # in form of cos theta
                # helices_angles.append(angle)
                # print(angle)
                out.write("{}, {}, {}\n".format(ii, count+1, angle))
            # helices_angles_all_frames.append(helices_angles)


def structure_prediction_run(protein):
    print(protein)
    protocol_list = ["awsemer", "frag", "er"]
    do = os.system
    cd = os.chdir
    cd(protein)
    # run = "frag"
    for protocol in protocol_list:
        do("rm -r " + protocol)
        do("mkdir -p " + protocol)
        do("cp -r {} {}/".format(protein, protocol))
        cd(protocol)
        cd(protein)
        # do("cp ~/opt/gremlin/protein/{}/gremlin/go_rnativeC* .".format(protein))
        do("cp ~/opt/gremlin/protein/{}/raptor/go_rnativeC* .".format(protein))
        fileName = protein + "_multi.in"
        backboneFile = "fix_backbone_coeff_" + protocol
        with fileinput.FileInput(fileName, inplace=True, backup='.bak') as file:
            for line in file:
                tmp = line.replace("fix_backbone_coeff_er", backboneFile)
                print(tmp, end='')
        cd("..")
        do("run.py -m 0 -n 20 {}".format(protein))
        cd("..")
    cd("..")
    # do("")


def check_and_correct_fragment_memory(fragFile="fragsLAMW.mem"):
    with open("tmp.mem", "w") as out:
        with open(fragFile, "r") as f:
            for i in range(4):
                line = next(f)
                out.write(line)
            for line in f:
                gro, _, i, n, _ = line.split()
                delete = False
                # print(gro, i, n)
                # name = gro.split("/")[-1]
                with open(gro, "r") as one:
                    next(one)
                    next(one)
                    all_residues = set()
                    for atom in one:
                        residue, *_ = atom.split()
                        # print(residue)
                        all_residues.add(int(residue))
                    for test in range(int(i), int(i)+int(n)):
                        if test not in all_residues:
                            print("ATTENTION", gro, i, n, "missing:",test)
                            delete = True
                if not delete:
                    out.write(line)
    os.system(f"mv {fragFile} fragsLAMW_back")
    os.system(f"mv tmp.mem {fragFile}")

def read_complete_temper_2(n=4, location=".", rerun=-1, qnqc=False, average_z=False, localQ=False, disReal=False, dis_h56=False, goEnergy=False, goEnergy3H=False, goEnergy4H=False):
    all_data_list = []
    for i in range(n):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file)
        lipid.columns = lipid.columns.str.strip()
        remove_columns = ['Steps']
        lipid = lipid.drop(remove_columns, axis=1)

        file = "rgs.{}.dat".format(i)
        rgs = pd.read_csv(location+file)
        rgs.columns = rgs.columns.str.strip()
        remove_columns = ['Steps']
        rgs = rgs.drop(remove_columns, axis=1)

        file = "energy.{}.dat".format(i)
        energy = pd.read_csv(location+file)
        energy.columns = energy.columns.str.strip()
        energy = energy[["AMH-Go", "Membrane", "Rg"]]
        file = "addforce.{}.dat".format(i)
        dis = pd.read_csv(location+file)
        dis.columns = dis.columns.str.strip()
        remove_columns = ['Steps', 'AddedForce', 'Dis12', 'Dis34', 'Dis56']
        dis.drop(remove_columns, axis=1,inplace=True)


        file = "wham.{}.dat".format(i)
        wham = pd.read_csv(location+file).assign(Run=i)
        wham.columns = wham.columns.str.strip()
        remove_columns = ['Rg', 'Tc']
        wham = wham.drop(remove_columns, axis=1)
        if qnqc:
            qc = pd.read_table(location+f"qc_{i}", names=["qc"])[1:].reset_index(drop=True)
            qn = pd.read_table(location+f"qn_{i}", names=["qn"])[1:].reset_index(drop=True)
            qc2 = pd.read_table(location+f"qc2_{i}", names=["qc2"])[1:].reset_index(drop=True)
            wham = pd.concat([wham, qn, qc, qc2],axis=1)
        # if average_z:
        #     z = pd.read_table(location+f"z_{i}.dat", names=["AverageZ"])[1:].reset_index(drop=True)
        #     wham = pd.concat([wham, z],axis=1)
        if disReal:
            tmp = pd.read_csv(location+f"distance_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            tmp.columns = tmp.columns.str.strip()
            wham = pd.concat([wham, tmp],axis=1)
        if dis_h56:
            tmp = pd.read_csv(location+f"distance_h56_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            tmp1 = pd.read_csv(location+f"distance_h12_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            tmp2 = pd.read_csv(location+f"distance_h34_{i}.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            tmp.columns = tmp.columns.str.strip()
            tmp1.columns = tmp1.columns.str.strip()
            tmp2.columns = tmp2.columns.str.strip()
            wham = pd.concat([wham, tmp, tmp1, tmp2],axis=1)
        if average_z:
            z = pd.read_csv(location+f"z_complete_{i}.dat")[1:].reset_index(drop=True)
            z.columns = z.columns.str.strip()
            wham = pd.concat([wham, z],axis=1)
        if localQ:
            all_localQ = pd.read_csv(location+f"localQ.{i}.csv")[1:].reset_index(drop=True)
            wham = pd.concat([wham, all_localQ], axis=1)
        if goEnergy:
            tmp = pd.read_csv(location+f"Go_{i}/goEnergy.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            tmp.columns = tmp.columns.str.strip()
            wham = pd.concat([wham, tmp],axis=1)
        if goEnergy3H:
            nEnergy = pd.read_csv(location+f"Go_3helix_{i}/goEnergy.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            nEnergy.columns = nEnergy.columns.str.strip()
            wham = pd.concat([wham, nEnergy],axis=1)
        if goEnergy4H:
            nEnergy = pd.read_csv(location+f"Go_4helix_{i}/goEnergy.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
            # print(tmp)
            nEnergy.columns = nEnergy.columns.str.strip()
            wham = pd.concat([wham, nEnergy],axis=1)
        data = pd.concat([wham, dis, energy, rgs, lipid], axis=1)

        # lipid = lipid[["Steps","Lipid","Run"]]
        all_data_list.append(data)
    data = pd.concat(all_data_list)
    file = f"../log{rerun}/log.lammps"
    temper = pd.read_table(location+file, skiprows=2, sep=' ')
    temper = temper.melt(id_vars=['Step'], value_vars=['T' + str(i) for i in range(n)], value_name="Temp", var_name="Run")
    temper["Run"] = temper["Run"].str[1:].astype(int)
    temper["Temp"] = "T" + temper["Temp"].astype(str)
#     print(temper)
#     print(wham)
    t2 = temper.merge(data, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
#     print(t2)
    t3 = t2.assign(TotalE=t2.Energy + t2.Lipid)
    return t3.sort_values(["Step", "Run"]).reset_index(drop=True)

def process_complete_temper_data_3(pre, data_folder, folder_list, rerun=-1, end=-1, n=12, bias="dis", qnqc=False, average_z=False, disReal=False, dis_h56=False, localQ=False, goEnergy=False, goEnergy3H=False, goEnergy4H=False, label=""):
    print("process temp data")
    dateAndTime = datetime.today().strftime('%d_%h_%H%M%S')
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        # this one only consider rerun >=0, for the case rerun=-1, move log.lammps to log0
        for i in range(rerun, end, -1):
            all_data_list = []
            for one_simulation in simulation_list:
                bias_num = one_simulation.split("_")[-1]
                print(bias_num, "!")

                location = one_simulation + f"/{i}/"
                print(location)
                data = read_complete_temper_2(location=location, n=n, rerun=i, qnqc=qnqc, average_z=average_z, localQ=localQ, disReal=disReal, dis_h56=dis_h56, goEnergy=goEnergy, goEnergy3H=goEnergy3H, goEnergy4H=goEnergy4H)
                print(data.shape)
                # remove_columns = ['Step', "Run"]
                # data = data.drop(remove_columns, axis=1)
                all_data_list.append(data.assign(BiasTo=bias_num))

            data = pd.concat(all_data_list).reset_index(drop=True)
            # if localQ:
            #     print("hi")
            # else:
            #     data.to_csv(os.path.join(pre, folder, f"data/rerun_{i}.csv"))
            # complete_data_list.append(data)
            #         temps = list(dic.keys())
            # complete_data = pd.concat(complete_data_list)
            name = f"rerun_{2*i}_{dateAndTime}.feather"
            data = data.reset_index(drop=True)
            data.query(f'Step > {2*i}e7 & Step <= {2*i+1}e7').reset_index(drop=True).to_feather(pre+folder+"/" + name)
            os.system("cp "+pre+folder+"/" + name + " "+data_folder+label+name)
            name = f"rerun_{2*i+1}_{dateAndTime}.feather"
            data = data.reset_index(drop=True)
            data.query(f'Step > {2*i+1}e7 & Step <= {2*i+2}e7').reset_index(drop=True).to_feather(pre+folder+"/" + name)
            os.system("cp "+pre+folder+"/" + name + " "+data_folder+label+name)



def move_data4(data_folder, freeEnergy_folder, folder_list, temp_dict_mode=1, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, biasName="dis", qnqc=False, average_z=0, chosen_mode=0):
    print("move data")
    # dic = {"T_defined":300, "T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
    if temp_dict_mode == 1:
        dic = {"T0":280, "T1":300, "T2":325, "T3":350, "T4":375, "T5":400, "T6":450, "T7":500, "T8":550, "T9":600, "T10":650, "T11":700}
    if temp_dict_mode == 2:
        dic = {"T0":280, "T1":290, "T2":300, "T3":315, "T4":335, "T5":355, "T6":380, "T7":410, "T8":440, "T9":470, "T10":500, "T11":530}
    if temp_dict_mode == 3:
        dic = {"T0":280, "T1":290, "T2":300, "T3":310, "T4":320, "T5":335, "T6":350, "T7":365, "T8":380, "T9":410, "T10":440, "T11":470}
    if temp_dict_mode == 4:
        dic = {"T0":300, "T1":335, "T2":373, "T3":417, "T4":465, "T5":519, "T6":579, "T7":645, "T8":720, "T9":803, "T10":896, "T11":1000}
    # read in complete.feather
    data_list = []
    for folder in folder_list:
        tmp = pd.read_feather(data_folder + folder +".feather")
        data_list.append(tmp)
    data = pd.concat(data_list)
    os.system("mkdir -p "+freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}")
    for bias, oneBias in data.groupby("BiasTo"):
        for tempSymbol, oneTempAndBias in oneBias.groupby("Temp"):
            temp = dic[tempSymbol]
            if float(temp) > 800:
                continue
            print(f"t_{temp}_{biasName}_{bias}.dat")
            if sample_range_mode == 0:
                queryCmd = 'Step > 0 & Step <= 1e7'
            if sample_range_mode == 1:
                queryCmd = 'Step > 1e7 & Step <= 2e7'
            elif sample_range_mode == 2:
                queryCmd ='Step > 2e7 & Step <= 3e7'
            elif sample_range_mode == 3:
                queryCmd ='Step > 3e7 & Step <= 4e7'
            elif sample_range_mode == 4:
                queryCmd ='Step > 4e7 & Step <= 5e7'
            elif sample_range_mode == 5:
                queryCmd ='Step > 5e7 & Step <= 6e7'
            elif sample_range_mode == 6:
                queryCmd ='Step > 6e7 & Step <= 7e7'
            elif sample_range_mode == 7:
                queryCmd ='Step > 7e7 & Step <= 8e7'
            elif sample_range_mode == -1:
                queryCmd ='Step > 4e7 & Step <= 6e7'
            if sample_range_mode == -2:
                tmp = oneTempAndBias.reset_index(drop=True)
            else:
                tmp = oneTempAndBias.query(queryCmd).reset_index()
            if average_z < 5:
                chosen_list = ["TotalE", "Qw", "Distance"]
            elif average_z == 5:
                chosen_list = ["TotalE", "Qw", "DisReal"]
                chosen_list += ["z_h6"]
            if average_z == 1:
                chosen_list += ["abs_z_average"]
            if average_z == 2 or average_z == 3:
                chosen_list += ["z_h6"]
            if average_z == 3:
                chosen_list += ["DisReal"]
            if average_z == 4:
                tmp["z_h5_and_h6"] = tmp["z_h5"] + tmp["z_h6"]
                chosen_list += ["z_h5_and_h6"]
                chosen_list += ["DisReal"]
            if average_z == 6:
                chosen_list = ["TotalE", "Qw", "DisReal"]
                tmp["z_h5_and_h6"] = tmp["z_h5"] + tmp["z_h6"]
                chosen_list += ["z_h5_and_h6"]
                chosen_list += ["z_h5"]
                chosen_list += ["z_h6"]
                chosen_list += ["Dis_h56"]
            if average_z == 7:
                chosen_list = ["TotalE", "Qw", "DisReal"]
                tmp["z_h56"] = tmp["z_h5"] + tmp["z_h6"]
                tmp["z_h14"] = tmp["z_h1"] + tmp["z_h2"] + tmp["z_h3"] + tmp["z_h4"]
                chosen_list += ["z_h14"]
                chosen_list += ["z_h56"]
                chosen_list += ["z_h5"]
                chosen_list += ["z_h6"]
                chosen_list += ["Dis_h12"]
                chosen_list += ["Dis_h34"]
                chosen_list += ["Dis_h56"]
            if chosen_mode == 0:
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg)
            if chosen_mode == 1:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
            if chosen_mode == 2:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg,
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg)
    #         print(tmp.count())
            if chosen_mode == 3:
                chosen_list += ["AMH-Go", "Lipid", "Membrane", "Rg"]
                chosen = tmp[chosen_list]
            if chosen_mode == 4:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
            if chosen_mode == 5:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_go_m=tmp.TotalE/10,
                                        TotalE_perturb_go_p=0,
                                        Go=tmp["AMH-Go"])
            if chosen_mode == 6:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH,
                                        TotalE_4=tmp.TotalE + tmp.AMH,
                                        TotalE_5=tmp.AMH)
            if chosen_mode == 7:
                chosen_list += ["Dis_h56"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_3H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_3H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_3H,
                                        TotalE_4=tmp.TotalE + tmp.AMH_3H,
                                        TotalE_5=tmp.TotalE + 0.1*tmp.AMH,
                                        TotalE_6=tmp.TotalE + 0.2*tmp.AMH)
            if chosen_mode == 8:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_4H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_4H,
                                        TotalE_4=tmp.TotalE + 0.1*tmp.AMH_3H,
                                        TotalE_5=tmp.TotalE + 0.2*tmp.AMH_3H,
                                        TotalE_6=tmp.TotalE + 0.5*tmp.AMH_3H)
            if chosen_mode == 9:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_4H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_4H)
                chosen = chosen.assign(TotalE_perturb_1go_m=chosen.TotalE_2 - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_1go_p=chosen.TotalE_2 + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_2lipid_m=chosen.TotalE_2 - tmp.Lipid,
                                        TotalE_perturb_2lipid_p=chosen.TotalE_2 + tmp.Lipid,
                                        TotalE_perturb_3mem_m=chosen.TotalE_2 - tmp.Membrane,
                                        TotalE_perturb_3mem_p=chosen.TotalE_2 + tmp.Membrane,
                                        TotalE_perturb_4rg_m=chosen.TotalE_2 - tmp.Rg,
                                        TotalE_perturb_4rg_p=chosen.TotalE_2 + tmp.Rg,
                                        TotalE_perturb_5go=tmp["AMH-Go"],
                                        TotalE_perturb_5lipid=tmp.Lipid,
                                        TotalE_perturb_5mem=tmp.Membrane,
                                        TotalE_perturb_5rg=tmp.Rg)
            if chosen_mode == 10:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 0.1*tmp.AMH_4H,
                                        TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_3=tmp.TotalE + 0.5*tmp.AMH_4H)
                chosen = chosen.assign(TotalE_perturb_1lipid_m1=chosen.TotalE_2 - 0.1*tmp.Lipid,
                                        TotalE_perturb_1lipid_p1=chosen.TotalE_2 + 0.1*tmp.Lipid,
                                        TotalE_perturb_2lipid_m2=chosen.TotalE_2 - 0.2*tmp.Lipid,
                                        TotalE_perturb_2lipid_p2=chosen.TotalE_2 + 0.2*tmp.Lipid,
                                        TotalE_perturb_3lipid_m3=chosen.TotalE_2 - 0.3*tmp.Lipid,
                                        TotalE_perturb_3lipid_p3=chosen.TotalE_2 + 0.3*tmp.Lipid,
                                        TotalE_perturb_4lipid_m4=chosen.TotalE_2 - 0.5*tmp.Lipid,
                                        TotalE_perturb_4lipid_p4=chosen.TotalE_2 + 0.5*tmp.Lipid,
                                        TotalE_perturb_5go=tmp["AMH-Go"],
                                        TotalE_perturb_5lipid=tmp.Lipid,
                                        TotalE_perturb_5mem=tmp.Membrane,
                                        TotalE_perturb_5rg=tmp.Rg)
            if chosen_mode == 11:
                # chosen_list += ["Dis_h56"]
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_1=tmp.TotalE + 1.1*0.1*tmp.AMH_4H + 0.1*tmp["AMH-Go"],
                                        TotalE_2=tmp.TotalE + 1.1*0.2*tmp.AMH_4H + 0.1*tmp["AMH-Go"],
                                        TotalE_3=tmp.TotalE + 1.1*0.5*tmp.AMH_4H + 0.1*tmp["AMH-Go"])
                chosen = chosen.assign(TotalE_perturb_1lipid_m1=chosen.TotalE_2 - 0.1*tmp.Lipid,
                                        TotalE_perturb_1lipid_p1=chosen.TotalE_2 + 0.1*tmp.Lipid,
                                        TotalE_perturb_2lipid_m2=chosen.TotalE_2 - 0.2*tmp.Lipid,
                                        TotalE_perturb_2lipid_p2=chosen.TotalE_2 + 0.2*tmp.Lipid,
                                        TotalE_perturb_3lipid_m3=chosen.TotalE_2 - 0.1*tmp.Membrane,
                                        TotalE_perturb_3lipid_p3=chosen.TotalE_2 + 0.1*tmp.Membrane,
                                        TotalE_perturb_4lipid_m4=chosen.TotalE_2 - 0.2*tmp.Membrane,
                                        TotalE_perturb_4lipid_p4=chosen.TotalE_2 + 0.2*tmp.Membrane,
                                        TotalE_perturb_5go=tmp["AMH-Go"],
                                        TotalE_perturb_5lipid=tmp.Lipid,
                                        TotalE_perturb_5mem=tmp.Membrane,
                                        TotalE_perturb_5rg=tmp.Rg)
            if chosen_mode == 12:
                chosen = tmp[chosen_list]
                # chosen["z_h56"] = (chosen["z_h5"] + chosen["z_h6"])/2
                chosen = chosen.assign(TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        z_h56=(tmp.z_h5 + tmp.z_h6)/2)
            if chosen_mode == 13:
                chosen_list += ["z_average"]
                chosen = tmp[chosen_list]
                # chosen["z_h56"] = (chosen["z_h5"] + chosen["z_h6"])/2
                force = 0.1
                chosen = chosen.assign(TotalE_2=tmp.TotalE + 0.2*tmp.AMH_4H - (tmp.DisReal - 25.1)*force,
                                        TotalE_3=tmp.TotalE - (tmp.DisReal - 25.1)*force,
                                        TotalE_4=tmp.TotalE + 0.2*tmp.AMH_4H,
                                        TotalE_5=tmp.TotalE + 0.2*tmp.AMH_4H - (tmp.DisReal)*force)
            chosen.to_csv(freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}/t_{temp}_{biasName}_{bias}.dat", sep=' ', index=False, header=False)

    # perturbation_table = {0:"original", 1:"m_go",
    #                       2:"p_go", 3:"m_lipid",
    #                       4:"p_lipid", 5:"m_mem",
    #                       6:"p_mem", 7:"m_rg", 8:"p_rg"}
def compute_average_z(dumpFile, outFile):
    # input dump, output z.dat
    z_list = []
    with open(outFile, "w") as f:
        a = read_lammps(dumpFile)
        for atoms in a:
            b = np.array(atoms)
            z = b.mean(axis=0)[2]
            z_list.append(z)
            f.write(str(z)+"\n")

def compute_average_z_2(dumpFile, outFile):
    # input dump, output z.dat

    helices_list = [(94,114), (147,168), (171, 192), (200, 217), (226, 241), (250, 269)]
    with open(outFile, "w") as f:
        a = read_lammps(dumpFile)
        f.write("z_average, abs_z_average, z_h1, z_h2, z_h3, z_h4, z_h5, z_h6\n")
        for atoms in a:
            b = np.array(atoms)
            z = b.mean(axis=0)[2]
            f.write(str(z)+ ", ")
            z = np.abs(b).mean(axis=0)[2]
            f.write(str(z)+ ", ")
            for count, (i,j) in enumerate(helices_list):
                i = i - 91
                j = j - 91
                z = np.mean(b[i:j], axis=0)[2]
                if count == 5:
                    f.write(str(z))
                else:
                    f.write(str(z)+ ", ")
            f.write("\n")

def read_simulation_2(location=".", i=-1, qnqc=False, average_z=False, localQ=False, disReal=False, **kwargs):
    file = "lipid.dat"
    lipid = pd.read_csv(location+file)
    lipid.columns = lipid.columns.str.strip()
    remove_columns = ['Steps']
    lipid = lipid.drop(remove_columns, axis=1)

    file = "rgs.dat"
    rgs = pd.read_csv(location+file)
    rgs.columns = rgs.columns.str.strip()
    remove_columns = ['Steps']
    rgs = rgs.drop(remove_columns, axis=1)

    file = "energy.dat"
    energy = pd.read_csv(location+file)
    energy.columns = energy.columns.str.strip()
    energy = energy[["AMH-Go", "Membrane", "Rg"]]
    file = "addforce.dat"
    dis = pd.read_csv(location+file)
    dis.columns = dis.columns.str.strip()
    remove_columns = ['Steps', 'AddedForce', 'Dis12', 'Dis34', 'Dis56']
    dis.drop(remove_columns, axis=1,inplace=True)


    file = "wham.dat"
    wham = pd.read_csv(location+file).assign(Run=i)
    wham.columns = wham.columns.str.strip()
    remove_columns = ['Rg', 'Tc']
    wham = wham.drop(remove_columns, axis=1)
    if qnqc:
        qc = pd.read_table(location+f"qc", names=["qc"])[1:].reset_index(drop=True)
        qn = pd.read_table(location+f"qn", names=["qn"])[1:].reset_index(drop=True)
        qc2 = pd.read_table(location+f"qc2", names=["qc2"])[1:].reset_index(drop=True)
        wham = pd.concat([wham, qn, qc, qc2],axis=1)
    # if average_z:
    #     z = pd.read_table(location+f"z_{i}.dat", names=["AverageZ"])[1:].reset_index(drop=True)
    #     wham = pd.concat([wham, z],axis=1)
    if disReal:
        tmp = pd.read_csv(location+f"distance.dat")[1:].reset_index(drop=True).drop('Steps', axis=1)
        # print(tmp)
        tmp.columns = tmp.columns.str.strip()
        wham = pd.concat([wham, tmp],axis=1)
    if average_z:
        z = pd.read_csv(location+f"z_complete.dat")[1:].reset_index(drop=True)
        z.columns = z.columns.str.strip()
        wham = pd.concat([wham, z],axis=1)
    if localQ:
        all_localQ = pd.read_csv(location+f"localQ.csv")[1:].reset_index(drop=True)
        wham = pd.concat([wham, all_localQ], axis=1)

    data = pd.concat([wham, dis, energy, rgs, lipid], axis=1)
    t3 = data.assign(TotalE=data.Energy + data.Lipid)
    return t3.reset_index(drop=True)

def read_folder(location, match="", **kwargs):
    runFolders = os.listdir(location+"/simulation")
    if match == "qbias":
        runFolders = [f for f in runFolders if re.match(r'qbias_[0-9]+', f)]
    else:
        runFolders = [f for f in runFolders if re.match(r'[0-9]+', f)]
    print(runFolders)
    data_list = []
    for run in runFolders:
        tmp = read_simulation_2(location+"/simulation/"+run+"/0/", **kwargs).assign(Run=run)
        data_list.append(tmp)
    return pd.concat(data_list).reset_index(drop=True)

def read_variable_folder(location, match="*_", **kwargs):
    variables = glob.glob(os.path.join(location, match))
    print(variables)
    data_list = []
    for variableFolder in variables:
        tmp = variableFolder.split("/")[-1]
        data_list.append(read_folder(variableFolder, **kwargs).assign(Folder=tmp))
    data = pd.concat(data_list)
    name = f"{datetime.today().strftime('%d_%h_%H%M%S')}.feather"
    data.reset_index(drop=True).to_feather(name)


def downloadPdb(pdb_list):
    os.system("mkdir -p original_pdbs")
    for pdb_id in pdb_list:
        pdb = f"{pdb_id.lower()[:4]}"
        pdbFile = pdb+".pdb"
        if not os.path.isfile("original_pdbs/"+pdbFile):
            pdbl = PDBList()
            name = pdbl.retrieve_pdb_file(pdb, pdir='.', file_format='pdb')
            os.system(f"mv {name} original_pdbs/{pdbFile}")



def cleanPdb(pdb_list, chain=None):
    os.system("mkdir -p cleaned_pdbs")
    for pdb_id in pdb_list:
        pdb = f"{pdb_id.lower()[:4]}"
        if chain is None:
            if len(pdb_id) == 5:
                Chosen_chain = pdb_id[4].upper()
            else:
                assert(len(pdb_id) == 4)
                Chosen_chain = "A"
        else:
            Chosen_chain = chain
        pdbFile = pdb+".pdb"
        # clean pdb
        fixer = PDBFixer(filename="original_pdbs/"+pdbFile)
        # remove unwanted chains
        chains = list(fixer.topology.chains())
        chains_to_remove = [i for i, x in enumerate(chains) if x.id not in Chosen_chain]
        fixer.removeChains(chains_to_remove)

        fixer.findMissingResidues()
        # add missing residues in the middle of a chain, not ones at the start or end of the chain.
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        # print(keys)
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        PDBFile.writeFile(fixer.topology, fixer.positions, open("cleaned_pdbs/"+pdbFile, 'w'))


def getAllChains(pdbFile):
    fixer = PDBFixer(filename=pdbFile)
    # remove unwanted chains
    chains = list(fixer.topology.chains())
    a = ""
    for i in chains:
        a += i.id
    return ''.join(sorted(set(a.upper().replace(" ", ""))))


def add_chain_to_pymol_pdb(location):
    # location = "/Users/weilu/Research/server/nov_2018/openMM/random_start/1r69.pdb"
    with open("tmp", "w") as out:
        with open(location, "r") as f:
            for line in f:
                info = list(line)
                if len(info) > 21:
                    info[21] = "A"
                out.write("".join(info))
    os.system(f"mv tmp {location}")















# ----------------------------depreciated---------------------------------------



def read_simulation(location):
    file = "lipid.dat"
    lipid = pd.read_csv(location+file)
    lipid.columns = lipid.columns.str.strip()

    file = "energy.dat"
    energy = pd.read_csv(location+file)
    energy.columns = energy.columns.str.strip()
    file = "addforce.dat"
    dis = pd.read_csv(location+file)
    dis.columns = dis.columns.str.strip()
#     remove_columns = ['AddedForce', 'Dis12', 'Dis34', 'Dis56']
    file = "rgs.dat"
    rgs = pd.read_csv(location+file)
    rgs.columns = rgs.columns.str.strip()
    file = "wham.dat"
    wham = pd.read_csv(location+file)
    wham.columns = wham.columns.str.strip()
    remove_columns = ['Rg', 'Tc']
    wham = wham.drop(remove_columns, axis=1)
    data = wham.merge(rgs, how='inner', left_on=["Steps"], right_on=["Steps"]).\
        merge(dis, how='inner', left_on=["Steps"], right_on=["Steps"]).\
        merge(energy, how='inner', left_on=["Steps"], right_on=["Steps"]).\
        merge(lipid, how='inner', left_on=["Steps"], right_on=["Steps"])
    data = data.assign(TotalE=data.Energy + data.Lipid)
    return data


def process_complete_temper_data_2(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis", qnqc=False, average_z=False, localQ=False):
    print("process temp data")
    dateAndTime = datetime.today().strftime('%d_%h_%H%M%S')
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        # this one only consider rerun >=0, for the case rerun=-1, move log.lammps to log0
        for i in range(rerun+1):
            all_data_list = []
            for one_simulation in simulation_list:
                bias_num = one_simulation.split("_")[-1]
                print(bias_num, "!")

                location = one_simulation + f"/{i}/"
                print(location)
                data = read_complete_temper_2(location=location, n=n, rerun=i, qnqc=qnqc, average_z=average_z, localQ=localQ)
                print(data.shape)
                # remove_columns = ['Step', "Run"]
                # data = data.drop(remove_columns, axis=1)
                all_data_list.append(data.assign(BiasTo=bias_num))

            data = pd.concat(all_data_list).reset_index(drop=True)
            # if localQ:
            #     print("hi")
            # else:
            #     data.to_csv(os.path.join(pre, folder, f"data/rerun_{i}.csv"))
            # complete_data_list.append(data)
            #         temps = list(dic.keys())
            # complete_data = pd.concat(complete_data_list)
            name = f"rerun_{i}_{dateAndTime}.feather"
            data.reset_index(drop=True).to_feather(pre+folder+"/" + name)
            os.system("cp "+pre+folder+"/" + name + " "+data_folder)

def move_data3(data_folder, freeEnergy_folder, folder, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, biasName="dis", qnqc=False, average_z=0, chosen_mode=0):
    print("move data")
    dic = {"T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
    # read in complete.feather
    data = pd.read_feather(data_folder + folder +".feather")
    os.system("mkdir -p "+freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}")
    for bias, oneBias in data.groupby("BiasTo"):
        for tempSymbol, oneTempAndBias in oneBias.groupby("Temp"):
            temp = dic[tempSymbol]
            if float(temp) > 800:
                continue
            print(f"t_{temp}_{biasName}_{bias}.dat")
            if sample_range_mode == 0:
                queryCmd = 'Step > 0 & Step <= 1e7'
            if sample_range_mode == 1:
                queryCmd = 'Step > 1e7 & Step <= 2e7'
            elif sample_range_mode == 2:
                queryCmd ='Step > 2e7 & Step <= 3e7'
            elif sample_range_mode == 3:
                queryCmd ='Step > 3e7 & Step <= 4e7'
            elif sample_range_mode == 4:
                queryCmd ='Step > 4e7 & Step <= 5e7'
            elif sample_range_mode == 5:
                queryCmd ='Step > 5e7 & Step <= 6e7'
            elif sample_range_mode == -1:
                queryCmd ='Step > 4e7 & Step <= 6e7'
            tmp = oneTempAndBias.query(queryCmd)
            chosen_list = ["TotalE", "Qw", "Distance"]
            if average_z == 1:
                chosen_list += ["abs_z_average"]
            if average_z == 2:
                chosen_list += ["z_h6"]
            if chosen_mode == 0:
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg)
            if chosen_mode == 1:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
    #         print(tmp.count())
            chosen.to_csv(freeEnergy_folder+"/"+sub_mode_name+f"/data_{sample_range_mode}/t_{temp}_{biasName}_{bias}.dat", sep=' ', index=False, header=False)


def move_data2(data_folder, freeEnergy_folder, folder, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, biasName="dis", qnqc=False, average_z=0, chosen_mode=0):
    print("move data")
    dic = {"T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
    # read in complete.feather
    data = pd.read_feather(data_folder + folder +".feather")
    os.system("mkdir -p "+freeEnergy_folder+folder+sub_mode_name+"/data")
    for bias, oneBias in data.groupby("BiasTo"):
        for tempSymbol, oneTempAndBias in oneBias.groupby("Temp"):
            temp = dic[tempSymbol]
            if float(temp) > 800:
                continue
            print(f"t_{temp}_{biasName}_{bias}.dat")
            if sample_range_mode == 0:
                queryCmd = 'Step > 1e7 & Step <= 2e7'
            elif sample_range_mode == 1:
                queryCmd ='Step > 2e7 & Step <= 3e7'
            elif sample_range_mode == 2:
                queryCmd ='Step > 3e7 & Step <= 4e7'
            elif sample_range_mode == 3:
                queryCmd ='Step > 4e7 & Step <= 5e7'
            elif sample_range_mode == 4:
                queryCmd ='Step > 5e7 & Step <= 6e7'
            elif sample_range_mode == -1:
                queryCmd ='Step > 4e7 & Step <= 6e7'
            tmp = oneTempAndBias.query(queryCmd)
            chosen_list = ["TotalE", "Qw", "Distance"]
            if average_z == 1:
                chosen_list += ["abs_z_average"]
            if average_z == 2:
                chosen_list += ["z_h6"]
            if chosen_mode == 0:
                chosen = tmp[chosen_list]
                chosen = chosen.assign(TotalE_perturb_mem_p=tmp.TotalE + kmem*tmp.Membrane,
                                        TotalE_perturb_mem_m=tmp.TotalE - kmem*tmp.Membrane,
                                        TotalE_perturb_lipid_p=tmp.TotalE + klipid*tmp.Lipid,
                                        TotalE_perturb_lipid_m=tmp.TotalE - klipid*tmp.Lipid,
                                        TotalE_perturb_go_p=tmp.TotalE + kgo*tmp["AMH-Go"],
                                        TotalE_perturb_go_m=tmp.TotalE - kgo*tmp["AMH-Go"],
                                        TotalE_perturb_rg_p=tmp.TotalE + krg*tmp.Rg,
                                        TotalE_perturb_rg_m=tmp.TotalE - krg*tmp.Rg)
            if chosen_mode == 1:
                chosen_list += ["Res" + str(i+1) for i in range(181)]
                chosen = tmp[chosen_list]
    #         print(tmp.count())
            chosen.to_csv(freeEnergy_folder+folder+sub_mode_name+f"/data/t_{temp}_{biasName}_{bias}.dat", sep=' ', index=False, header=False)
    # chosen

def make_metadata_2(cwd=".", k=1000.0, temps_list=["450"]):
    files = glob.glob("../../data/*")
    kconstant = k
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            bias = tmp.split("_")[3]
            # print(tmp)
            # if int(float(dis)) > 150:
            #     continue
            if t in temps_list:
                target = "{} {} {} {}\n".format(oneFile, t, kconstant, bias)
                out.write(target)

def make_metadata(k=1000.0, temps_list=["450"]):
    cwd = os.getcwd()
    files = glob.glob("../data/*")
    kconstant = k
    with open("metadatafile", "w") as out:
        for oneFile in files:
            tmp = oneFile.split("/")[-1].replace('.dat', '')
            t = tmp.split("_")[1]
            bias = tmp.split("_")[3]
            # print(tmp)
            # if int(float(dis)) > 150:
            #     continue
            if t in temps_list:
                target = "../{} {} {} {}\n".format(oneFile, t, kconstant, bias)
                out.write(target)
def read_complete_temper(n=4, location=".", rerun=-1, qnqc=False, average_z=False, localQ=False):
    all_lipid_list = []
    for i in range(n):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file).assign(Run=i)
        lipid.columns = lipid.columns.str.strip()
        # lipid = lipid[["Steps","Lipid","Run"]]
        all_lipid_list.append(lipid)
    lipid = pd.concat(all_lipid_list)

    all_rgs_list = []
    for i in range(n):
        file = "rgs.{}.dat".format(i)
        rgs = pd.read_csv(location+file).assign(Run=i)
        rgs.columns = rgs.columns.str.strip()
        # lipid = lipid[["Steps","Lipid","Run"]]
        all_rgs_list.append(rgs)
    rgs = pd.concat(all_rgs_list)

    all_energy_list = []
    for i in range(n):
        file = "energy.{}.dat".format(i)
        energy = pd.read_csv(location+file).assign(Run=i)
        energy.columns = energy.columns.str.strip()
        energy = energy[["Steps", "AMH-Go", "Membrane", "Rg", "Run"]]
        all_energy_list.append(energy)
    energy = pd.concat(all_energy_list)

    all_dis_list = []
    for i in range(n):
        file = "addforce.{}.dat".format(i)
        dis = pd.read_csv(location+file).assign(Run=i)
        dis.columns = dis.columns.str.strip()
        remove_columns = ['AddedForce', 'Dis12', 'Dis34', 'Dis56']
        dis.drop(remove_columns, axis=1,inplace=True)
        all_dis_list.append(dis)
    dis = pd.concat(all_dis_list)

    all_wham_list = []
    for i in range(n):
        file = "wham.{}.dat".format(i)
        wham = pd.read_csv(location+file).assign(Run=i)
        wham.columns = wham.columns.str.strip()
        remove_columns = ['Rg', 'Tc']
        wham = wham.drop(remove_columns, axis=1)
        if qnqc:
            qc = pd.read_table(location+f"qc_{i}", names=["qc"])[1:].reset_index(drop=True)
            qn = pd.read_table(location+f"qn_{i}", names=["qn"])[1:].reset_index(drop=True)
            qc2 = pd.read_table(location+f"qc2_{i}", names=["qc2"])[1:].reset_index(drop=True)
            wham = pd.concat([wham, qn, qc, qc2],axis=1)
        if average_z:
            z = pd.read_table(location+f"z_{i}.dat", names=["AverageZ"])[1:].reset_index(drop=True)
            wham = pd.concat([wham, z],axis=1)
        if localQ:
            all_localQ = pd.read_csv(location+f"localQ.{i}.csv")[1:].reset_index(drop=True)
            wham = pd.concat([wham, all_localQ], axis=1)
        all_wham_list.append(wham)
    wham = pd.concat(all_wham_list)
    if rerun == -1:
        file = "../log.lammps"
    else:
        file = f"../log{rerun}/log.lammps"
    temper = pd.read_table(location+file, skiprows=2, sep=' ')
    temper = temper.melt(id_vars=['Step'], value_vars=['T' + str(i) for i in range(n)], value_name="Temp", var_name="Run")
    temper["Run"] = temper["Run"].str[1:].astype(int)
    temper["Temp"] = "T" + temper["Temp"].astype(str)
    t2 = temper.merge(wham, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t3 = t2.merge(dis, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t4 = t3.merge(lipid, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t5 = t4.merge(energy, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t6 = t5.merge(rgs, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t6 = t6.assign(TotalE=t6.Energy + t6.Lipid)
    return t6
def process_complete_temper_data(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis", qnqc=False, average_z=False, localQ=False):
    print("process temp data")
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        complete_data_list = []
        for one_simulation in simulation_list:
            bias_num = one_simulation.split("_")[-1]
            print(bias_num, "!")
            all_data_list = []
            if rerun == -1:
                location = one_simulation + "/0/"
                print(location)
                data = read_complete_temper(location=location, n=n, rerun=rerun, qnqc=qnqc, average_z=average_z, localQ=localQ)
                # remove_columns = ['Step', "Run"]
                # data = data.drop(remove_columns, axis=1)
                all_data_list.append(data)
            else:
                for i in range(rerun+1):
                    location = one_simulation + f"/{i}/"
                    print(location)
                    data = read_complete_temper(location=location, n=n, rerun=i, qnqc=qnqc, average_z=average_z, localQ=localQ)
                    # remove_columns = ['Step', "Run"]
                    # data = data.drop(remove_columns, axis=1)
                    all_data_list.append(data)

            data = pd.concat(all_data_list).assign(BiasTo=bias_num).reset_index(drop=True)
            if localQ:
                print("hi")
            else:
                data.to_csv(os.path.join(pre, folder, f"data/bias_num.csv"))
            complete_data_list.append(data)
    #         temps = list(dic.keys())
        complete_data = pd.concat(complete_data_list)
        name = f"{datetime.today().strftime('%d_%h_%H%M%S')}.feather"
        complete_data.reset_index(drop=True).to_feather(pre+folder+"/" + name)
        os.system("cp "+pre+folder+"/" + name + " "+data_folder+folder+".feather")
def read_temper(n=4, location=".", rerun=-1, qnqc=False):
    all_lipid_list = []
    for i in range(n):
        file = "lipid.{}.dat".format(i)
        lipid = pd.read_csv(location+file).assign(Run=i)
        lipid.columns = lipid.columns.str.strip()
        lipid = lipid[["Steps","Lipid","Run"]]
        all_lipid_list.append(lipid)
    lipid = pd.concat(all_lipid_list)

    all_energy_list = []
    for i in range(n):
        file = "energy.{}.dat".format(i)
        energy = pd.read_csv(location+file).assign(Run=i)
        energy.columns = energy.columns.str.strip()
        energy = energy[["Steps", "AMH-Go", "Membrane", "Rg", "Run"]]
        all_energy_list.append(energy)
    energy = pd.concat(all_energy_list)

    all_dis_list = []
    for i in range(n):
        file = "addforce.{}.dat".format(i)
        dis = pd.read_csv(location+file).assign(Run=i)
        dis.columns = dis.columns.str.strip()
        remove_columns = ['AddedForce', 'Dis12', 'Dis34', 'Dis56']
        dis.drop(remove_columns, axis=1,inplace=True)
        all_dis_list.append(dis)
    dis = pd.concat(all_dis_list)

    all_wham_list = []
    for i in range(n):
        file = "wham.{}.dat".format(i)
        wham = pd.read_csv(location+file).assign(Run=i)
        wham.columns = wham.columns.str.strip()
        remove_columns = ['Rg', 'Tc']
        wham = wham.drop(remove_columns, axis=1)
        if qnqc:
            qc = pd.read_table(location+f"qc_{i}", names=["qc"])[1:].reset_index(drop=True)
            qn = pd.read_table(location+f"qn_{i}", names=["qn"])[1:].reset_index(drop=True)
            qc2 = pd.read_table(location+f"qc2_{i}", names=["qc2"])[1:].reset_index(drop=True)
            wham = pd.concat([wham,qn, qc, qc2],axis=1)
        all_wham_list.append(wham)
    wham = pd.concat(all_wham_list)
    if rerun == -1:
        file = "../log.lammps"
    else:
        file = f"../log{rerun}/log.lammps"
    temper = pd.read_table(location+file, skiprows=2, sep=' ')
    temper = temper.melt(id_vars=['Step'], value_vars=['T' + str(i) for i in range(n)], value_name="Temp", var_name="Run")
    temper["Run"] = temper["Run"].str[1:].astype(int)
    temper["Temp"] = "T" + temper["Temp"].astype(str)
    t2 = temper.merge(wham, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t3 = t2.merge(dis, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t4 = t3.merge(lipid, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t5 = t4.merge(energy, how='inner', left_on=["Step", "Run"], right_on=["Steps", "Run"]).sort_values('Step').drop('Steps', axis=1)
    t6 = t5.assign(TotalE=t5.Energy + t5.Lipid)
    return t6


def process_temper_data(pre, data_folder, folder_list, rerun=-1, n=12, bias="dis", qnqc=False):
    print("process temp data")
    for folder in folder_list:
        simulation_list = glob.glob(pre+folder+f"/simulation/{bias}_*")
        print(pre+folder+f"/simulation/{bias}_*")
        os.system("mkdir -p " + pre+folder+"/data")
        for one_simulation in simulation_list:
            bias_num = one_simulation.split("_")[-1]
            print(bias_num, "!")
            if rerun == -1:
                location = one_simulation + "/0/"
                try:
                    data = read_temper(location=location, n=n, qnqc=qnqc)
                    # remove_columns = ['Step', "Run"]
                    # data = data.drop(remove_columns, axis=1)
                    data.reset_index().to_feather(pre+folder+"/data/"+f"{bias}{bias_num}.feather")
                except:
                    print("Unexpected error:", sys.exc_info()[0])
                    print("notrun?", dis)
            else:
                all_data_list = []
                for i in range(rerun):
                    location = one_simulation + f"/{i}/"
                    try:
                        data = read_temper(location=location, n=n, rerun=i, qnqc=qnqc)
                        # remove_columns = ['Step', "Run"]
                        # data = data.drop(remove_columns, axis=1)
                        all_data_list.append(data)
                    except:
                        print("Unexpected error:", sys.exc_info()[0])
                        print("notrun?", bias_num)
                try:
                    data = pd.concat(all_data_list)
                    data.reset_index(drop=True).to_feather(pre+folder+"/data/"+f"{bias}{bias_num}.feather")
                except:
                    print("Unexpected error:", sys.exc_info()[0])
                    print("not data?", bias_num)
    #         temps = list(dic.keys())
        os.system("mv "+pre+folder+"/data "+data_folder+folder)

def move_data(data_folder, freeEnergy_folder, folder, sub_mode_name="", kmem=0.2, klipid=0.1, kgo=0.1, krg=0.2, sample_range_mode=0, bias="dis"):
    print("move data")
    os.system("mkdir -p "+freeEnergy_folder+folder+sub_mode_name+"/data")
    dis_list = glob.glob(data_folder+folder+f"/{bias}*.feather")
    for dis_file in dis_list:
        dis = dis_file.split("/")[-1].replace(bias, '').replace('.feather', '')
        print(dis)
        t6 = pd.read_feather(dis_file)
        remove_columns = ['index']
        t6 = t6.drop(remove_columns, axis=1)
        t6 = t6.assign(TotalE_perturb_mem_p=t6.TotalE + kmem*t6.Membrane)
        t6 = t6.assign(TotalE_perturb_mem_m=t6.TotalE - kmem*t6.Membrane)
        t6 = t6.assign(TotalE_perturb_lipid_p=t6.TotalE + klipid*t6.Lipid)
        t6 = t6.assign(TotalE_perturb_lipid_m=t6.TotalE - klipid*t6.Lipid)
        t6 = t6.assign(TotalE_perturb_go_p=t6.TotalE + kgo*t6["AMH-Go"])
        t6 = t6.assign(TotalE_perturb_go_m=t6.TotalE - kgo*t6["AMH-Go"])
        t6 = t6.assign(TotalE_perturb_rg_p=t6.TotalE + krg*t6.Rg)
        t6 = t6.assign(TotalE_perturb_rg_m=t6.TotalE - krg*t6.Rg)
        # t6["TotalE"] = 0
        dic = {"T0":350, "T1":400, "T2":450, "T3":500, "T4":550, "T5":600, "T6":650, "T7":700, "T8":750, "T9":800, "T10":900, "T11":1000}
        temps = list(dic.values())

        def convert(x):
            return dic[x]
        t6["Temp"] = t6["Temp"].apply(convert)

        for temp in temps:
            if temp > 800:
                continue
            if sample_range_mode == 0:
                tmp = t6.query('Temp=="{}"& Step > 1e7 & Step <= 2e7'.format(temp))
            elif sample_range_mode == 1:
                tmp = t6.query('Temp=="{}"& Step > 2e7 & Step <= 3e7'.format(temp))
            elif sample_range_mode == 2:
                tmp = t6.query('Temp=="{}"& Step > 3e7 & Step <= 4e7'.format(temp))
            tmp.to_csv(freeEnergy_folder+folder+sub_mode_name+"/data/t_{}_{}_{}.dat".format(temp, bias, dis), sep=' ', index=False, header=False)
# def pick_out_and_show():
#     protein_list = ["1occ", "1pv6", "2bl2", "2bg9", "1j4n", "1py6"]
#     for protein in protein_list:
#         frames = [pd.read_csv("{}/awsemer/simulation/{}/0/wham.dat".format(protein, run)).assign(Run=run) for run in range(20)]
#         result = pd.concat(frames)
#         answer = result.iloc[result[' Qw'].argsort()].iloc[-1]
#         print(protein, answer.Steps, answer.Run)
#         os.chdir("{}/awsemer/simulation/{}/0/".format(protein, int(answer.Run)))
#         os.system("show.py --frame {} {} -p".format(int(answer.Steps/4000), protein))
#         os.chdir("../../../../../")
