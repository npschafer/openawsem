#!/usr/bin/env python3
import pandas as pd
import mdtraj as md
import numpy as np
import time
import openmm
import openawsem
import pytest
import functools
from pathlib import Path

"""
compares OpenAWSEM implementation of "old" hbond energy from LAMMPS AWSEM-MD commit cea754f 
(https://github.com/adavtyan/awsemmd/tree/cea754f1208fde6332d4d0f1cae3212bf7e8afbb) to energies computed by the LAMMPS code
"""


PROTEINS = ["6rb9","2ohx_A","8j47_IGECA"]#["1le5","2onv_6chains","2y3j_6chains","3loz_6chains","3nve_6chains","3nhc_6chains","3ow9_6chains","4r0p_6chains","2lnq","2l8x","7umq"]
COLUMNS = ["Beta","Pap",]# "Helical"]
PLATFORMS = ['Reference', 'CPU', 'OpenCL', 'CUDA']
data_path = Path('tests')/'data'/'test_implementation_of_lammps_hbond_energies'


def single_run(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        return result, elapsed_time*1000
    return wrapper

def repeated_run(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        total_time = 0
        executions = 0
        # Run the function repeatedly until the total time is at least 10 seconds
        while total_time < 20:
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            elapsed_time = end_time - start_time
            total_time += elapsed_time
            executions += 1
        # Calculate average time (in milliseconds) per execution
        average_time = (total_time / executions) * 1000
        return result, average_time
    return wrapper


@single_run
def time_once(func):
    return func()

@repeated_run
def time_many(func):
    return func()


def set_up_forces(oa, protein, force_name=None):
    #Define all forces using lambda to delay execution of the setup.
    all_forces = {
        "RamaSSWeight": lambda: openawsem.functionTerms.basicTerms.rama_ssweight_term(oa, k_rama_ssweight=2*8.368, ssweight_file=data_path/f'{protein}-ssweight'),
        "Beta1": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_1_old(oa,k_beta=0.5*4.184,ssweight_file=data_path/f'{protein}-ssweight'), 
        "Beta2": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_2_old(oa,k_beta=0.5*4.184,ssweight_file=data_path/f'{protein}-ssweight'), 
        "Beta3": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_3_old(oa,k_beta=0.5*4.184,ssweight_file=data_path/f'{protein}-ssweight'), 
        "Helical": lambda: openawsem.functionTerms.hydrogenBondTerms.helical_term(oa,forceGroup=29),
        "Pap": lambda: openawsem.functionTerms.hydrogenBondTerms.pap_term_old(oa,k_pap=0.5*4.184,ssweight_file=data_path/f'{protein}-ssweight',forceGroup=28),
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

    forceGroupTable = {"Beta": {23,24,25}, "Pap": {28}, "Helical": {29},} 

    termEnergies = pd.DataFrame(columns=["Step"] + COLUMNS)

    for step in range(len(pdb_trajectory)):
        simulation.context.setPositions(pdb_trajectory.openmm_positions(step))
        e = []
        for term in COLUMNS:
            g = forceGroupTable[term]
            state = simulation.context.getState(getEnergy=True, groups=g)
            # note that we use kcal/mol instead of kJ/mol because we're comparing to lammps, which uses kcal/mol
            termEnergy = state.getPotentialEnergy().value_in_unit(openawsem.unit_definitions.kilocalories_per_mole)
            e.append(termEnergy)
        termEnergies.loc[step] = [step] + e

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
        tolerance = 1e-2 # note the higher tolerance than the 1e-5 used in test_energies.py
        for protein in PROTEINS:
            calculated_energies = analyzed_data(protein, platform)
            saved_energies = pd.read_csv(data_path/f'{protein}_energies.csv')

            assert column in calculated_energies.columns, f"Column {column} not found in calculated energies for protein {protein} on platform {platform}"
            assert column in saved_energies.columns, f"Column {column} not found in saved energies for protein {protein} on platform {platform}"
            assert np.allclose(calculated_energies[column], saved_energies[column], atol=tolerance), f"Energy terms comparison failed for protein {protein} on column {column} on platform {platform}"


if __name__ == '__main__':
    # not executed by pytest
    pass
