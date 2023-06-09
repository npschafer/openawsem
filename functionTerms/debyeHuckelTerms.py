try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
import numpy as np

def debye_huckel_term(self, k_dh=4.15*4.184, forceGroup=30, screening_length=1.0, chargeFile=None):
        # screening_length (in the unit of nanometers)
        print("Debye Huckel term is ON")
        k_dh *= self.k_awsem*0.1
        k_screening = 1.0
        # screening_length = 1.0  # (in the unit of nanometers)
        min_seq_sep = 1
        dh = CustomBondForce(f"{k_dh}*charge_i*charge_j/r*exp(-{k_screening}*r/{screening_length})")
        dh.addPerBondParameter("charge_i")
        dh.addPerBondParameter("charge_j")
        structure_interactions_dh = []
        if chargeFile is None:
            for i in range(self.nres):
                for j in range(i+min_seq_sep,self.nres):
                    charge_i = 0.0
                    charge_j = 0.0
                    if self.seq[i] == "R" or self.seq[i]=="K":
                        charge_i = 1.0
                    if self.seq[i] == "D" or self.seq[i]=="E":
                        charge_i = -1.0
                    if self.seq[j] == "R" or self.seq[j]=="K":
                        charge_j = 1.0
                    if self.seq[j] == "D" or self.seq[j]=="E":
                        charge_j = -1.0
                    if charge_i*charge_j!=0.0:
                        cb_atom_i = self.cb[i]
                        if cb_atom_i == -1:
                            cb_atom_i = self.ca[i]  # if mutated, and CB isn't found, then use CA instead
                        cb_atom_j = self.cb[j]
                        if cb_atom_j == -1:
                            cb_atom_j = self.ca[j]  # if mutated, and CB isn't found, then use CA instead
                        structure_interactions_dh.append([cb_atom_i, cb_atom_j, [charge_i, charge_j]])
                        # print([self.seq[i], self.seq[j],self.cb[i], self.cb[j], [charge_i, charge_j]])
        else:
            chargeInfo = np.loadtxt(chargeFile, dtype=[('index', int), ('charge', float)])
            for i in range(self.nres):
                charge_i = chargeInfo[i][1]
                for j in range(i+min_seq_sep,self.nres):
                    charge_j = chargeInfo[j][1]
                    if charge_i*charge_j!=0.0:
                        cb_atom_i = self.cb[i]
                        if cb_atom_i == -1:
                            cb_atom_i = self.ca[i]  # if mutated, and CB isn't found, then use CA instead
                        cb_atom_j = self.cb[j]
                        if cb_atom_j == -1:
                            cb_atom_j = self.ca[j]  # if mutated, and CB isn't found, then use CA instead
                        structure_interactions_dh.append([cb_atom_i, cb_atom_j, [charge_i, charge_j]])
        for structure_interaction_dh in structure_interactions_dh:
            dh.addBond(*structure_interaction_dh)
        dh.setForceGroup(forceGroup)
        return dh
