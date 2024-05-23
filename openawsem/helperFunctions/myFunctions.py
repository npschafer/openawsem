#!python
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import subprocess
import glob
import re
import copy

from openawsem.helperFunctions.myFunctions_helper import *
import numpy as np
import pandas as pd
import fileinput
from itertools import product
from Bio.PDB.PDBParser import PDBParser

from Bio.PDB import PDBList
from pdbfixer import PDBFixer
try:
    from openmm.app import PDBFile
except ModuleNotFoundError:
    from simtk.openmm.app import PDBFile

# compute cross Q for every pdb pair in one folder
# parser = argparse.ArgumentParser(description="Compute cross q")
# parser.add_argument("-m", "--mode",
#                     type=int, default=1)

# args = parser.parse_args()

def getFromTerminal(CMD):
    return subprocess.Popen(CMD,stdout=subprocess.PIPE,shell=True).communicate()[0].decode()


def expand_grid(dictionary):
    return pd.DataFrame([row for row in product(*dictionary.values())],
                        columns=dictionary.keys())


def duplicate_pdb(From, To, offset_x=0, offset_y=0, offset_z=0, new_chain="B"):
    with open(To, "w") as out:
        with open(From, "r") as f:
            for line in f:
                tmp = list(line)
                if len(tmp) < 26:
                    out.write(line)
                    continue
                atom = line[0:4]
                atomSerialNumber = line[6:11]
                atomName = line[12:16]
                atomResidueName = line[17:20]
                chain = line[21]
                residueNumber = line[22:26]
                # change chain A to B
                # new_chain = "B"

                if atom == "ATOM":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # add 40 to the x
                    new_x = x + offset_x
                    new_y = y + offset_y
                    new_z = z + offset_z

                    tmp[21] = new_chain
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
                    all_residues = []
                    for atom in one:
                        residue, resType, atomType, *_ = atom.split()
                        # print(residue, resType, atomType)
                        if atomType == "CA":
                            all_residues.append(int(residue))
                    all_residues = np.array(all_residues)
                    for test in range(int(i), int(i)+int(n)):
                        if (test == all_residues).sum() > 1:
                            # In rare case, one res id may have two different possible residues.
                            # on example, pdb 3vpg. chain A, res id 220.
                            # ATOM   1467  N   ARG A 220A      9.151 -20.984  46.737  1.00 31.30           N
                            # ATOM   1468  CA  ARG A 220A      9.120 -19.710  46.027  1.00 31.52           C
                            # ATOM   1469  C   ARG A 220A      9.768 -19.832  44.650  1.00 33.58           C
                            # ATOM   1470  O   ARG A 220A     10.552 -18.973  44.240  1.00 28.91           O
                            # ATOM   1471  CB  ARG A 220A      9.853 -18.641  46.847  1.00 31.58           C
                            # ATOM   1472  CG  ARG A 220A      9.181 -18.295  48.168  1.00 33.55           C
                            # ATOM   1473  CD  ARG A 220A      7.834 -17.651  47.916  1.00 34.70           C
                            # ATOM   1474  NE  ARG A 220A      7.959 -16.526  46.994  1.00 43.05           N
                            # ATOM   1475  CZ  ARG A 220A      6.931 -15.906  46.425  1.00 46.69           C
                            # ATOM   1476  NH1 ARG A 220A      5.691 -16.300  46.683  1.00 39.12           N
                            # ATOM   1477  NH2 ARG A 220A      7.144 -14.898  45.590  1.00 41.15           N
                            # ATOM   1478  N   ALA A 220B      9.429 -20.901  43.936  1.00 33.78           N
                            # ATOM   1479  CA  ALA A 220B      9.979 -21.153  42.608  1.00 32.13           C
                            # ATOM   1480  C   ALA A 220B      9.944 -19.933  41.692  1.00 30.71           C
                            # ATOM   1481  O   ALA A 220B      9.050 -19.088  41.787  1.00 28.56           O
                            # ATOM   1482  CB  ALA A 220B      9.234 -22.310  41.951  1.00 35.20           C
                            print("ATTENTION", gro, i, n, "duplicate:",test)
                            delete = True
                        if test not in all_residues:
                            print("ATTENTION", gro, i, n, "missing:",test)
                            delete = True
                if not delete:
                    out.write(line)
    os.system(f"mv {fragFile} fragsLAMW_back")
    os.system(f"mv tmp.mem {fragFile}")

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


def downloadPdb(pdb_list, membrane_protein=False, location="original_pdbs/"):
    print("Download from server")
    os.system(f"mkdir -p {location}")
    for pdb_id in pdb_list:
        pdb = f"{pdb_id.lower()[:4]}"
        pdbFile = pdb+".pdb"
        fileLocation = os.path.join(location, pdbFile)
        if not os.path.isfile(fileLocation):
            if membrane_protein:
                # os.system(f"wget http://pdbtm.enzim.hu/data/database/fn/{pdbFile}.gz")
                # os.system(f"gunzip {pdbFile}.gz")
                os.system(f"wget https://opm-assets.storage.googleapis.com/pdb/{pdbFile}")
                os.system(f"mv {pdbFile} {fileLocation}")
            else:
                pdbl = PDBList()
                name = pdbl.retrieve_pdb_file(pdb, pdir='.', file_format='pdb')
                os.system(f"mv {name} {fileLocation}")
            os.system("rm -r obsolete")



def cleanPdb(pdb_list, chain=None, source=None, toFolder="cleaned_pdbs", formatName=False, 
                removeDNAchains=True, verbose=False, removeTwoEndsMissingResidues=True, addMissingResidues=True, removeHeterogens=True, keepIds=False):
    os.system(f"mkdir -p {toFolder}")
    for pdb_id in pdb_list:
        # print(chain)
        print(pdb_id)
        # pdb = f"{pdb_id.lower()[:4]}"
        # pdbFile = pdb+".pdb"
        if formatName:
            pdb = f"{pdb_id.lower()[:4]}"
        else:
            pdb = pdb_id
        pdbFile = pdb + ".pdb"
        if source is None:
            fromFile = os.path.join("original_pdbs", pdbFile)
        elif source[-4:] == ".pdb":
            fromFile = source
        else:
            fromFile = os.path.join(source, pdbFile)

        if verbose:
            print('Fixing PDB using PDBFixer')
            print(os.getcwd())
            print(fromFile)
            
        # clean pdb
        fixer = PDBFixer(filename=fromFile)

        try:
            fixer = PDBFixer(filename=fromFile)
        except Exception as inst:
            print(inst)
            print(f"{fromFile} not found. skipped")
            continue
        
        if verbose:
            print('Removing unwanted chains')
        
        # remove unwanted chains
        chains = list(fixer.topology.chains())
        print(chains)
        if chain is None:  # 'None' means deafult is chain A unless specified.
            if len(pdb_id) >= 5:
                Chosen_chain = pdb_id[4]
                # Chosen_chain = pdb_id[4].upper()
            else:
                assert(len(pdb_id) == 4)
                Chosen_chain = "A"
        elif chain == "-1" or chain == -1:
            Chosen_chain = getAllChains(fromFile, removeDNAchains=removeDNAchains)
            print(f"Chains: {Chosen_chain}")
        elif chain == "first":
            Chosen_chain = chains[0].id
        else:
            Chosen_chain = chain
            
        chains_to_remove = [i for i, x in enumerate(chains) if x.id not in Chosen_chain]
        fixer.removeChains(chains_to_remove)

        if verbose:
            print('Adding Missing residues')
        
        fixer.findMissingResidues()
        # add missing residues in the middle of a chain, not ones at the start or end of the chain.
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        if verbose:
            print("chains to remove", chains_to_remove)
            print("missing residues: ", keys)
        if not addMissingResidues:
            for key in list(keys):
                del fixer.missingResidues[key]
        else:
            if removeTwoEndsMissingResidues:
                for key in list(keys):
                    chain_tmp = chains[key[0]]
                    if key[1] == 0 or key[1] == len(list(chain_tmp.residues())):
                        del fixer.missingResidues[key]

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        if removeHeterogens:
            fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        try:
            fixer.addMissingAtoms()
        except:
            print("Unable to add missing atoms")
            continue
        fixer.addMissingHydrogens(7.0)
        PDBFile.writeFile(fixer.topology, fixer.positions, open(os.path.join(toFolder, pdbFile), 'w'), keepIds=keepIds)


def getAllChains(pdbFile, removeDNAchains=True):
    fixer = PDBFixer(filename=str(pdbFile))
    # we only want pdb chains, ligands or DNA chain will be ignored here.
    fixer.removeHeterogens(keepWater=False)
    # remove unwanted chains
    chains = list(fixer.topology.chains())
    a = ""

    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    rnaResidues = ['A', 'G', 'C', 'U', 'I']
    dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']
    for c in chains:
        if removeDNAchains and np.alltrue([a.name in dnaResidues for a in c.residues()]):
            print(f"chain {c.id} is a DNA chain. it will be removed")
            continue
        if c.id in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789':
            a += c.id
    # return ''.join(sorted(set(a.upper().replace(" ", ""))))
    return ''.join(sorted(set(a.replace(" ", ""))))

def add_chain_to_pymol_pdb(location):
    # location = "/Users/weilu/Research/server/nov_2018/openMM/random_start/1r69.pdb"
    with open("tmp", "w") as out:
        with open(location, "r") as f:
            lines = f.readlines()
            if len(lines) < 1:
                print("your extended.pdb is empty. please check. (a possible cause is that the pymol is not installed)")
            for line in lines:
                info = list(line)
                if len(info) > 21:
                    if line[:4] == "ATOM":
                        info[21] = "A"
                out.write("".join(info))
    os.system(f"mv tmp {location}")


def get_seq_dic(fasta="../crystal_structure.fasta"):
    seq_dic = {}
    chain = None
    with open(fasta) as f:
        for line in f:
            if line[0] == ">":
                assert line[:19] == ">CRYSTAL_STRUCTURE:"
                if chain is not None:
                    seq_dic[chain] = seq
                chain = line[19]
                seq = ""
            else:
                seq += line.replace("\n", "")
        seq_dic[chain] = seq
    return seq_dic


def seq_length_from_pdb(fileLocation, chains):
    data = []
    parser = PDBParser()
    structure = parser.get_structure('X', fileLocation)
    chain_start_residue_index = 1
    for c in structure.get_chains():
        chain_name = c.get_id()
        if chain_name in chains:
            seq_len = len(list(c.get_residues()))
            print(chain_name, seq_len)
            data.append((chain_name, chain_start_residue_index, seq_len))
            chain_start_residue_index += seq_len
    return data


def get_frame(file="movie.pdb", to="last_frame.pdb", frame=-1):
    # default is last frame.
    # if you want first, please set frame to 1.
    a = open(file).read().split("ENDMDL")
    assert a[-1] == "\nEND\n"
    with open(to, "w") as out:
        out.write(a[frame-1])



def convert_openMM_to_standard_pdb(fileName="last_frame.pdb", seq_dic=None, back=True):
    code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
            "ARG" : "R", "LYS" : "K", "MET" : "M", "CYS" : "C",
            "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
            "TRP" : "W", "ASP" : "D", "GLU" : "E", "ASN" : "N",
            "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V"}
    inv_code_map = {v: k for k, v in code.items()}
    if seq_dic is None:
        seq_dic = get_seq_dic()
    if back and os.path.exists(fileName+".bak"):
        backup = '.bak'
    else:
        backup = ''
    with fileinput.FileInput(fileName, inplace=True, backup=backup) as file:
        for line in file:
            if len(line) >= 4 and line[:4] == "END\n":
                continue
            if len(line) > 25:
                if line[:6] == "REMARK":
                    continue
                i = int(line[22:26])
                chain = line[21]

                tmp = list(line)
                if "".join(tmp[17:20]) in ["IGL", "NGP", "IPR"]:
                    res = seq_dic[chain][i-1]
                    tmp[17:20] = inv_code_map[res]
                    if line[:6] == "HETATM":
                        tmp[:6] = "ATOM  "
                else:
                    # no change
                    pass
                print("".join(tmp), end='')
            else:
                print(line, end='')


# def relocate(location):
#     # location = "/Users/weilu/Research/server/april_2019/iterative_optimization_new_set_with_frag/all_simulations/1fc2/1fc2"
#     fileLocation = location + "/frags.mem"
#     # pre = location + "/../"
#     pre = location
#     os.system(f"mkdir -p {pre}/fraglib")
#     a = pd.read_csv(fileLocation, skiprows=4, sep=" ", names=["location", "i", "j", "sep", "w"])
#     b = a["location"].unique()
#     for l in b:
#         out = os.system(f"cp {l} {pre}/fraglib/")
#         if out != 0:
#             print(f"!!Problem!!, {l}")

def relocate(fileLocation="frags.mem", toLocation="fraglib"):
    # location = "/Users/weilu/Research/server/april_2019/iterative_optimization_new_set_with_frag/all_simulations/1fc2/1fc2"
    # fileLocation = location + "/frags.mem"
    # toLocation
    print(os.getcwd())
    os.system(f"mkdir -p {toLocation}")
    a = pd.read_csv(fileLocation, skiprows=4, sep=" ", names=["location", "i", "j", "sep", "w"])
    b = a["location"].unique()
    for l in b:
        cmd = f"cp {l} {toLocation}/"
        out = subprocess.Popen(cmd, shell=True).wait()
        if out != 0:
            print(f"!!Problem!!, {l}")

def replace(TARGET, FROM, TO):
    os.system("sed -i.bak 's@{}@{}@g' {}".format(FROM,TO,TARGET))




def get_PDB_length(pdbFileLocation):
    from Bio.PDB.PDBParser import PDBParser
    # pdbFileLocation = '/Users/weilu/Research/database/chosen/T0869-D1.pdb'
    structure = PDBParser().get_structure("a", pdbFileLocation)
    return len(list(structure.get_residues()))

def pdbToFasta(pdb, pdbLocation, fastaFile, chains="A"):
    import textwrap
    from Bio.PDB.PDBParser import PDBParser
    three_to_one = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
                    'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I',
                    'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P',
                    'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

    # pdb = "1r69"
    # pdbLocation = "/Users/weilu/Research/server/may_2019/family_fold/1r69.pdb"
    # chains = "A"
    # fastaFile = "/Users/weilu/Research/server/may_2019/family_fold/1r69.fasta"
    s = PDBParser().get_structure("X", pdbLocation)
    m = s[0]  # model 0
    seq = ""
    with open(fastaFile, "w") as out:
        for chain in chains:
            out.write(f">{pdb.upper()}:{chain.upper()}\n")
            c = m[chain]
            chain_seq = ""
            for residue in c:
                is_regular_res = residue.has_id('CA') and residue.has_id('O')
                res_id = residue.get_id()[0]
                if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
                    residue_name = residue.get_resname()
                    chain_seq += three_to_one[residue_name]
            out.write("\n".join(textwrap.wrap(chain_seq, width=80))+"\n")
            seq += chain_seq
    return seq


def read_hydrophobicity_scale(seq, tableLocation, isNew=False):
    seq_dataFrame = pd.DataFrame({"oneLetterCode":list(seq)})
    # HFscales = pd.read_table("~/opt/small_script/Whole_residue_HFscales.txt")
    # print(f"reading hydrophobicity scale table from {tableLocation}/Whole_residue_HFscales.txt")
    HFscales = pd.read_csv(f"{tableLocation}/Whole_residue_HFscales.txt", sep="\t")
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

def create_zim(fastaFile, tableLocation, isNew=True):
    # print("creating zim file for membrane potential")
    seq = ""
    with open(fastaFile, "r") as f:
        for line in f:
            if line[0] == ">":
                pass
            else:
                # print(line)
                seq += line.strip()
    data = read_hydrophobicity_scale(seq, tableLocation, isNew=isNew)
    z = data["DGwoct"].values
    np.savetxt("zim", z, fmt="%.2f")

def read_fasta(fastaFile):
    seq = ""
    with open(fastaFile, "r") as f:
        for line in f:
            if line[0] == ">":
                pass
            else:
                # print(line)
                seq += line.strip()
    return seq

def create_extended_pdb_from_fasta(filename, output_file_name="output.pdb"):

    # Standard distances and angles for the backbone atoms
    # proline is handled later
    CA_CA =3.8
    CA_C = 1.52/np.sqrt(3/2)
    CA_N = 1.45/np.sqrt(3/2)
    C_O = 1.23
    CA_CB = 1.53/np.sqrt(3/2)
    C_N = 1.33/np.sqrt(3/2)
    t=1/np.sqrt(2)

    # list of lists, where each sublist contains all coordinates (N, CA, C, O, CB)
    # of a single residue
    coords = []

    # Load the sequence from the FASTA file
    # Skip the name line and concatenate sequence lines
    with open(filename, "r") as file:
        lines = file.readlines()
    sequence = ''.join(line.strip() for line in lines[1:])
    print(sequence)

    # build extended chain residue by residue
    # start at origin and build in positive x direction until hitting proline.
    # place the proline appropriately, 
    # then translate all placed residues such that
    # the initial building process along the positive x axis produces a valid structure
    #
    # parameters for placing residues are inspired by pymol's peptide building functionality,
    # but are less sophisticated than pymol and will not reproduce the pymol structure
    #
    # start by defining the coordinates to use when placing the first residue at origin
    # these may need to be changed at proline resets because the backbone oxygens
    # have to alter orientation, so this will depend on whether an even or odd number
    # of residues was placed since the previous reset
    default_coord = np.array([[-1*CA_N, 0, -t*CA_N],  #N
                              [0, 0, 0],              #CA
                              [1*CA_C,0,-t*CA_C],     #C
                              [1*CA_C,0,-t*CA_C-C_O], #O
                              [0,-1*CA_CB,t*CA_CB]])  #CB
    default_coord[:,2]+=t*(CA_C+CA_N)/2
    # initialize before first iteration
    coord = np.array([default_coord[:,0]-CA_CA/np.sqrt(3**2+t**2)*3,
                     -default_coord[:,1],
                     -default_coord[:,2]]).T
    # the main loop
    for resid, residue in enumerate(sequence):
        # build and adjust chain as needed
        if sequence[resid] != "P":
            # standard build procedure
            coord[:,0]+=CA_CA/np.sqrt(3**2+t**2)*3
            coord[:,1]=-coord[:,1]
            coord[:,2]=-coord[:,2]
            # update chain with new residue
            coords.append(copy.deepcopy(coord))
        elif sequence[resid] == "P":
            # proline build procedure
            # as modification of standard build procedure
            coord[:,0]+=CA_CA/np.sqrt(3**2+t**2)*3 - 0.21
            coord[:,1]=-coord[:,1] + 0.55*((-1)**(resid+1)) # flip if 1st, 3rd, etc
            coord[:,2]=-coord[:,2]
            # now we need to rotate CB and O about C(previous residue)->N(proline) axis
            # it's easier to code this as an approximate translation rather than a real rotation
            coord[4,2] += 2.5*((-1)**(resid+1))
            coord[3,1] -= 2.23*((-1)**(resid+1))
            coord[2,1] -= 1.23*((-1)**(resid+1))
            # update chain with new residue
            coords.append(copy.deepcopy(coord))
            # reset chain
            #
            # fix the x coordinate by:
            # 1. subtracting what we added for proline
            for index,coord in enumerate(coords):
                coord[:,0] -= CA_CA/np.sqrt(3**2+t**2)*3 - 0.21
                coords[index] = coord
            # 2. subtracting the non-proline amount for every residue up to the previous proline
            for resid2 in range((len(sequence)-resid)+1,(len(sequence)+1)):
                #print(f'sequence: {sequence}')
                #print(f'-resid2: {-resid2}')
                if sequence[-resid2] == "P":
                    break
                else:
                    for index,coord in enumerate(coords):
                        coord[:,0] -= CA_CA/np.sqrt(3**2+t**2)*3
                        coords[index] = coord
            # fix the y and z coordinates by translating everything such that proline C
            # aligns with C from default_coord
            translation = default_coord[2,:] - coords[-1][2,:]
            for index,coord in enumerate(coords):
                coord += translation
                coords[index] = coord
            # rotate around x axis to adjust backbone alignment
            # it is important to do this before rotating in xy plane
            # but after translating y and z coordinates!
            #
            #first, temporarily translate such that y,z coordinates of proline C are 0
            old_y = copy.deepcopy([coord[2,1] for coord in coords])
            old_z = copy.deepcopy([coord[2,2] for coord in coords]) 
            for index,coord in enumerate(coords):
                #coord[:,1] -= coord[2,1] 
                #coord[:,2] -= coord[2,2] 
                coord[:,1] -= old_y[-1]
                coord[:,2] -= old_z[-1]
                coords[index] = coord
            # now do the rotation
            yz_theta = -39.1*np.pi/180 #0#21.6#39.4#-110#-90#-50#-30#-70.6
            yz_rotation_matrix = np.array([[np.cos(yz_theta),-np.sin(yz_theta)],[np.sin(yz_theta),np.cos(yz_theta)]])
            for index,coord in enumerate(coords):
                for atom_index in range(coord.shape[0]):
                    coord[atom_index,1:] = np.dot(yz_rotation_matrix,coord[atom_index,1:].T).T
                coords[index] = coord
            # now put it back where we found it
            for index,coord in enumerate(coords):
                #coord[:,1] += old_y[index]
                #coord[:,2] += old_z[index]
                coord[:,1] += old_y[-1]
                coord[:,2] += old_z[-1]
                coords[index] = coord
            # Now we need to rotate about the proline C-O bond axis
            #
            # first, translate such that proline C is at origin
            old = copy.deepcopy([coord[2,:] for coord in coords])
            for index,coord in enumerate(coords):
                coord -= coords[-1][2,:]
                coords[index] = coord
            # get C-O axis, which is just the oxygen coordinates
            CO_axis = coords[-1][3,:].T / np.linalg.norm(coords[-1][3,:].T,ord=2)
            # define rotation matrix and perform rotation
            CO_theta = 45*np.pi/180
            CO_rotation_matrix = np.array([[np.cos(CO_theta) + (CO_axis[0]**2)*(1-np.cos(CO_theta)),
                                                CO_axis[0]*CO_axis[1]*(1-np.cos(CO_theta)) - CO_axis[2]*np.sin(CO_theta),
                                                CO_axis[0]*CO_axis[2]*(1-np.cos(CO_theta)) + CO_axis[1]*np.sin(CO_theta)],
                                           [CO_axis[1]*CO_axis[0]*(1-np.cos(CO_theta)) + CO_axis[2]*np.sin(CO_theta),
                                                np.cos(CO_theta) + (CO_axis[1]**2)*(1-np.cos(CO_theta)),
                                                CO_axis[1]*CO_axis[2]*(1-np.cos(CO_theta)) - CO_axis[0]*np.sin(CO_theta)],
                                           [CO_axis[2]*CO_axis[0]*(1-np.cos(CO_theta)) - CO_axis[1]*np.sin(CO_theta),
                                                CO_axis[2]*CO_axis[1]*(1-np.cos(CO_theta)) + CO_axis[0]*np.sin(CO_theta),
                                                np.cos(CO_theta) + (CO_axis[2]**2)*(1-np.cos(CO_theta))]])
            for index,coord in enumerate(coords):
                for atom_index in range(coord.shape[0]):
                    coord[atom_index,:] = np.dot(CO_rotation_matrix,coord[atom_index,:].T).T
                coords[index] = coord
            # now put it back where we found it
            for index, coord in enumerate(coords):
                coord += old[-1]
                coords[index] = coord
            # after all of this, we have proline overlapping with the next residue after proline
            # so we'll translate 1 residue in x direction again
            for index, coord in enumerate(coords):
                coord[:,0] -= CA_CA/np.sqrt(3**2+t**2)*3
                coords[index] = coord
            # reset coord before continuing build
            if resid % 2 == 1: # the next one will be even
                coord = np.array([default_coord[:,0]-CA_CA/np.sqrt(3**2+t**2)*3,
                         -default_coord[:,1],
                         -default_coord[:,2]]).T
            elif resid % 2 == 0: # the next one will be odd
                coord = np.array([default_coord[:,0]-CA_CA/np.sqrt(3**2+t**2)*3,
                          default_coord[:,1],
                          default_coord[:,2]]).T
   
    # write coordinates to pdb file
    one_to_three={'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE',
                  'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 
                  'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN', 'R':'ARG',
                  'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}
    with open(output_file_name,'w') as f:
        for i in range(len(coords)): # each element of coords is the coordinates of a single residue
            coord = coords[i]
            for j, name in enumerate(['N','CA','C','O','CB']):
                # atom numbering is weird but this is fixed later
                f.write(f'ATOM  {j*5+1:>5} {name:^4} {one_to_three[sequence[i]]:<3} {"A"}{i:>4}    {coord[j][0]:>8.3f}{coord[j][1]:>8.3f}{coord[j][2]:>8.3f}\n')
        f.write("END\n")



















    '''
    proline_mask = []
    for one_letter_code in sequence:
        if one_letter_code == "P":
            proline_mask.append(1)
        else:
            proline_mask.append(0)
    with open(output_file_name, "w") as f:
        # For each residue in the sequence, calculate the position of the backbone atoms
        for resid, residue in enumerate(sequence):
            for i,name in enumerate(['N','CA','C','O','CB']):
                pdb_line = f'ATOM  {i*5+1:>5} {name:^4} {one_to_three[residue]:<3} {"A"}{resid:>4}    {coord[i][0]:>8.3f}{coord[i][1]:>8.3f}{coord[i][2]:>8.3f}\n' 
                f.write(pdb_line)
            # f.write("ATOM  %5d  N   %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+1, residue, i+1, coord[0][0], coord[0][1], coord[0][2]))
            # f.write("ATOM  %5d  CA  %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+2, residue, i+1, coord[1][0], coord[1][1], coord[1][2]))
            # f.write("ATOM  %5d  C   %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+3, residue, i+1, coord[2][0], coord[2][1], coord[2][2]))
            # f.write("ATOM  %5d  O   %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+4, residue, i+1, coord[3][0], coord[3][1], coord[3][2]))
            # f.write("ATOM  %5d  CB  %s A %3d    %8.3f%8.3f%8.3f\n" % (i*5+5, residue, i+1, coord[4][0], coord[4][1], coord[4][2]))
            
            # Calculate new coordinates for the next residue
            # this is unnecessary if we're on the last residue
            if resid < len(sequence)-1:
                # start with local modifications that affect construction of the residue,
                # in addition to global translation of the residue in the x dimension
                if sequence[resid] != "P" and sequence[resid+1] != "P":
                    coord[:,0]+=CA_CA/np.sqrt(3**2+t**2)*3
                    coord[:,1]=-coord[:,1]
                    coord[:,2]=-coord[:,2]
                elif sequence[resid] != "P" and sequence[resid+1] == "P":
                    coord[:,0]+=CA_CA/np.sqrt(3**2+t**2)*3 - 0.21
                    coord[:,1]=-coord[:,1] + 0.55
                    coord[:,2]=-coord[:,2]
                    # now we need to rotate CB and O about C(previous residue)->N(proline) axis
                    # it's easier to code this as an approximate translation rather than a real rotation
                    coord[4,2] += 2.5
                    coord[3,1] -= 2.23     
                    coord[2,1] -= 1.23
                elif sequence[resid] == "P" and sequence[resid+1] != "P": 
                    # we want to use non-proline parameters as a starting point, so we'll
                    # revert to older coord
                    ############################################################################
                    coord[:,0]+=CA_CA/np.sqrt(3**2+t**2)*3 -0.62 
                    coord[:,1]=-coord[:,1] -0.34
                    coord[:,2]=-coord[:,2] -1.09
                    #coord[2:4,0] -= .5
                    #coord[2:4,1] -= .5
                    #coord[2:4,2] -= .5
                    #############################################################################
                elif sequnece[resid] == "P" and sequence[resid+1] == "P":
                    raise ValueError("Two or more consecutive proline residues not supported")
                else:
                    raise AssertionError
                # rotate coord depending on the number of prolines that have been encountered
                # so far, including resid+1
                # this is necessary because prolines change the angle of the backbone
                num_prolines = sum(proline_mask[:resid+1]) # up to and including current resid
                                                           # because this iteration sets the coordinates for resid+1
                rotation_matrix = np.linalg.matrix_power(np.array([[np.cos(-46.8*np.pi/180),-np.sin(-46.8*np.pi/180)],
                        [np.sin(-46.8*np.pi/180),np.cos(-46.8*np.pi/180)]]), num_prolines)
                # this rotation matrix is appropriate for left-multiplication of a column vector [x,y]
                for atom_index in range(1,5):
                    # we'll fix the nitrogen atom as the origin for our rotation
                    temp = coord[atom_index,:2].reshape([2,1]) - coord[0,:2].reshape([2,1])
                    coord[atom_index,:2] = np.dot(rotation_matrix,temp).T + coord[0,:2]
                
                
                
                
                
                
                # for simplicity, we'll assume prolines and non-prolines have equal length   
                      # and their backbones are oriented in the same direction, except for when
                      # a proline follows a non-proline
                    

                        # global translation in the y and z dimensions, plus corrections to x dimension
                # to account for angled growth of chain after P
                # and combined effects of more than one P
        f.write("END\n")
        '''
