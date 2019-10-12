from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np

def debye_huckel_term(self, k_dh=4.15*4.184, forceGroup=30, screening_length = 1.0):
        # screening_length (in the unit of nanometers)
        print("Debye Huckel term is ON")
        k_dh *= self.k_awsem*0.1
        k_screening = 1.0


        dh = CustomBondForce(f"{k_dh}*charge_i*charge_j/r*exp(-{k_screening}*r/{screening_length})")
        dh.addPerBondParameter("charge_i")
        dh.addPerBondParameter("charge_j")
        structure_interactions_dh = []
        for i in range(self.nres):
            for j in range(i+1,self.nres):
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
                    structure_interactions_dh.append([self.cb[i], self.cb[j], [charge_i, charge_j]])
                    # print([self.seq[i], self.seq[j],self.cb[i], self.cb[j], [charge_i, charge_j]])
        for structure_interaction_dh in structure_interactions_dh:
            dh.addBond(*structure_interaction_dh)
        dh.setForceGroup(forceGroup)
        return dh
