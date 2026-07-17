#!/usr/bin/env python3

'''
While test_energies.py compares to reference energies that were
computed by openawsem at a particular point in time, this 
script compares to reference energies computed by lammps. 
Since the default precision of the coordinates written to
dump.lammpstrj is too low to get precise energy agreement for
some strong and sharply varying terms, the precision used by
lammps must be increased and different test functions must be 
used to load the trajectory and pass coordinates to openawsem
for energy evaluation. 
'''

import pandas as pd
import mdtraj as md
import numpy as np
import time
import openmm
import openawsem
import pytest
import functools
from pathlib import Path

PROTEINS = ["1brs",] #"1mbn", "1ubq", "2lyz", "2lzm"] # could add these others at a later date
COLUMNS = ["Con", "Chain", "Chi", "Excluded", "Rama",] #"Contact", "Fragment", "Membrane", "ER", "TBM_Q", "Beta", "Pap", "Helical"] # could add these others at a later date
PLATFORMS = ['Reference', 'CPU', 'OpenCL', 'CUDA']
data_path = Path('tests')/'data'

def set_up_forces(oa, protein, force_name=None):
    #Define all forces using lambda to delay execution of the setup.
    all_forces = {
        "Con": lambda: openawsem.functionTerms.basicTerms.con_term(oa, forceGroup=4),        
        "Rama": lambda: openawsem.functionTerms.basicTerms.rama_term(oa),
        "Contact": lambda: openawsem.functionTerms.contactTerms.contact_term(oa),
        "Chain": lambda: openawsem.functionTerms.basicTerms.chain_term(oa, forceGroup=5),
        "Chi": lambda: openawsem.functionTerms.basicTerms.chi_term(oa, forceGroup=6),
        "Excluded": lambda: openawsem.functionTerms.basicTerms.excl_term(oa, forceGroup=7),
        "RamaProline": lambda: openawsem.functionTerms.basicTerms.rama_proline_term(oa),
        "RamaSSWeight": lambda: openawsem.functionTerms.basicTerms.rama_ssweight_term(oa, k_rama_ssweight=2*8.368, ssweight_file=data_path/f'{protein}-ssweight'),
        #"Beta1": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_1(oa,ssweight_file=data_path/f'{protein}-ssweight'),
        #"Beta2": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_2(oa,ssweight_file=data_path/f'{protein}-ssweight'),
        #"Beta3": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_3(oa,ssweight_file=data_path/f'{protein}-ssweight'),
        #"Helical": lambda: openawsem.functionTerms.hydrogenBondTerms.helical_term(oa),
        #"Pap1": lambda: openawsem.functionTerms.hydrogenBondTerms.pap_term_1(oa,ssweight_file=data_path/f'{protein}-ssweight'),
        #"Pap2": lambda: openawsem.functionTerms.hydrogenBondTerms.pap_term_2(oa,ssweight_file=data_path/f'{protein}-ssweight'),
        #"FragmentMemory": lambda: openawsem.functionTerms.templateTerms.fragment_memory_term(oa, frag_file_list_file=data_path/f'{protein}-single_frags.mem', npy_frag_table=data_path/f'{protein}-single_frags.npy', UseSavedFragTable=False),
        #"DebyeHuckel": lambda: openawsem.functionTerms.debyeHuckelTerms.debye_huckel_term(oa, chargeFile=data_path/f'{protein}-charge.txt'),
    }
    forces = []
    if force_name:
        if force_name not in all_forces:
            raise ValueError(f"Force {force_name} is not recognized.")
        force = all_forces[force_name]()
        forces.append(force)
    else:
        for force_name, force_func in all_forces.items():
            force = force_func()
            forces.append(force)

    return forces

def analyze(protein, simulation_platform):
    # load sequence and get atom indices in lammps topologies corresponding 
    #     to starts and ends of different chains
    seq = []
    with open(data_path/f"{protein}-crystal_structure.fasta",'r') as f:
        current = ''
        for line in f:
            if line[0] != '>':
                current += line.strip()
            else:
                if current != '':
                    seq.append(current)
                    current = ''
    seq.append(current)
    print(f'fasta parser found the following chains: {seq}')
    concat_seq = ''.join(seq)
    n_atoms_per_chain = np.array([3*len(chain_seq) for chain_seq in seq]) # all lammps awsemmd residues have exactly 3 particles
    n_atoms = np.sum(n_atoms_per_chain)
    chain_starts = [1]
    for increment in n_atoms_per_chain[:-1]:
        chain_starts.append(chain_starts[-1]+increment)
    chain_starts = np.array(chain_starts)
    chain_ends = [n_atoms_per_chain[0]]
    for increment in n_atoms_per_chain[1:]:
        chain_ends.append(chain_ends[-1]+increment)
    chain_ends = np.array(chain_ends)
    # load high-precision lammps coordinates into a correctly shaped array
    #     lammps only logs CA, CB, and O, but we have to add positions for the virtual sites
    coordses = []
    coords = None
    with open(data_path/f'{protein}_lammps/dump.lammpstrj','r') as f: # high-precision coordinates
        for counter, line in enumerate(f):
            if counter % (n_atoms+9) < 9: # header lines
                if coords:
                    coordses.append(coords)
                coords = []
            else:
                info = line.strip().split(" ")
                try:
                    atom_id = int(info[0])
                except:
                    raise
                atom_type = int(info[1])
                xyz = [float(number) for number in info[2:]]
                if atom_type == 5: # the "HB" atom type from awsemmd is never included in openawsem topologies
                    continue
                elif np.any((0<=atom_id-chain_starts) & (atom_id-chain_starts<=2)): # first residue, which has a different virtual site arrangement than interior ones
                    if atom_type == 1: 
                        coords.append(xyz) # these are our CA coords
                    elif atom_type == 3:
                        coords.append([0,0,0]) # these are our C' virtual site coords, which will be fixed later
                        coords.append(xyz) # these are our O coords
                    elif atom_type == 4:
                        coords.append(xyz) # these are our CB coords
                    else:
                        raise AssertionError()
                elif np.any((chain_ends-2<=atom_id) & (atom_id<=chain_ends)): # last residue, which has a different virtual site arrangement than interior ones
                    if atom_type == 1:
                        coords.append([0,0,0]) # these are our N virtual site coords, which will be fixed later
                        if concat_seq[(atom_id-1)//3] != "P": 
                            coords.append([0,0,0]) # H virtual site coords, which only exist if not proline
                        coords.append(xyz) # CA coords of the last residue
                    elif atom_type == 3:
                        coords.append(xyz) # O coords of the last residue
                    elif atom_type == 4:
                        coords.append(xyz) # CB coords of the last residue
                    else:
                        raise AssertionError()
                else: # we are on an interior residue
                    if atom_type == 1:
                        coords.append([0,0,0]) # these are our N virtual site coords, which will be fixed later
                        if concat_seq[(atom_id-1)//3] != "P":
                            coords.append([0,0,0]) # H virtual site coords, which only exist if not proline
                        coords.append(xyz) # CA coords 
                    elif atom_type == 3:
                        coords.append([0,0,0]) # these are our C' virtual site coords, which will be fixed later
                        coords.append(xyz) # O coords
                    elif atom_type == 4:
                        coords.append(xyz) # CB coords
                    else:
                        raise AssertionError('unexpected else block')
    coordses.append(coords)
    coordses = np.array(coordses)/10 # don't forget the unit conversion!
    # normal analysis setup
    chain = openawsem.helperFunctions.myFunctions.getAllChains(data_path/f"{protein}-crystal_structure.pdb")
    seq = openawsem.helperFunctions.myFunctions.read_fasta(data_path/f"{protein}-crystal_structure.fasta")
    pdb_trajectory = md.load(data_path/f'{protein}-movie.dcd', top=data_path/f"{protein}-openmmawsem.pdb")
    oa = openawsem.OpenMMAWSEMSystem(data_path/f"{protein}-openmmawsem.pdb",
                                     chains=chain,
                                     k_awsem=1.0,
                                     xml_filename=openawsem.xml,
                                     seqFromPdb=seq,
                                     includeLigands=False)
    forces = set_up_forces(oa, protein)
    oa.addForcesWithDefaultForceGroup(forces)
    platform = openmm.Platform.getPlatformByName(simulation_platform)
    integrator = openmm.LangevinIntegrator(300*openawsem.unit_definitions.kelvin, 1/openawsem.unit_definitions.picosecond, 2*openawsem.unit_definitions.femtoseconds)
    simulation = openmm.app.Simulation(oa.pdb.topology, oa.system, integrator, platform)
    forceGroupTable = {"Con": 4, "Chain":5, "Chi":6, "Excluded":7, 
        "Rama": 21, "Contact": 22, "Fragment": 23, "Membrane": 24, "ER": 25, "TBM_Q": 26, "Beta": 27, "Pap": 28, "Helical": 29,
        "Q": 1, "Rg": 2, "Qc": 3, "Helix_orientation": 18, "Pulling": 19}
    # use high-precision coordinates when setting positions
    termEnergies = pd.DataFrame(columns=["Step"] + COLUMNS)
    for frame_index in range(coordses.shape[0]):
        simulation.context.setPositions(coordses[frame_index,:,:])
        simulation.context.computeVirtualSites() 
        e = []
        for term in COLUMNS:
            g = {forceGroupTable[term]} if forceGroupTable[term] != -1 else -1
            state = simulation.context.getState(getEnergy=True, groups=g)
            termEnergy = state.getPotentialEnergy().value_in_unit(openawsem.unit_definitions.kilojoule_per_mole if term in COLUMNS else openawsem.unit_definitions.kilocalories_per_mole)
            e.append(termEnergy)
        termEnergies.loc[frame_index] = [frame_index] + e

    return termEnergies

# Cache to store analyzed data
analyzed_data_cache = {}

@pytest.fixture(scope="session")
def analyzed_data():
    def get_data(protein, simulation_platform):
        # Create a unique key for each protein and platform combination
        key = (protein, simulation_platform)

        # If the data has already been computed, return it from the cache
        if key in analyzed_data_cache:
            return analyzed_data_cache[key]

        # Otherwise, compute the data and store it in the cache
        analyzed_data_cache[key] = analyze(protein, simulation_platform)
        return analyzed_data_cache[key]

    # Return the function that accesses data, either from cache or by computing
    return get_data

@pytest.mark.parametrize("column", COLUMNS)
@pytest.mark.parametrize("platform", PLATFORMS)
class TestEnergyTerms:
    def test_energy_term(self, platform, column, analyzed_data):
        tolerance = 5e-6 
        # Right now, the tolerance is limited by the fact that we only saved
        # 5 decimal places for some of the lammps energies. We could probably
        # make the tolerance lower if we had lammps write energies in double-
        # precision.
        for protein in PROTEINS:
            calculated_energies = analyzed_data(protein, platform)
            # lammps uses kcal/mol, but our openawsem energies will be calculated in kJ/mol,
            # so we convert here
            saved_energies = 4.184*pd.read_csv(data_path/f'{protein}_lammps/{protein}-lammps_energies.csv')

            assert column in calculated_energies.columns, f"Column {column} not found in calculated energies for protein {protein} on platform {platform}"
            assert column in saved_energies.columns, f"Column {column} not found in saved energies for protein {protein} on platform {platform}"
            assert np.allclose(calculated_energies[column], saved_energies[column], atol=tolerance), f"Energy terms comparison failed for protein {protein} on column {column} on platform {platform}"


if __name__ == '__main__':
    pass # in the main test script, test_energies.py, we offer an efficiency benchmark, 
         # but this one we want to keep simple for clarity