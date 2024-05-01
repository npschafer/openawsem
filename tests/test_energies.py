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

PROTEINS = ["1brs", "1mbn", "1ubq", "2lyz", "2lzm"]
COLUMNS = ["Backbone", "Rama", "Contact", "Fragment", "Membrane", "ER", "TBM_Q", "Beta", "Pap", "Helical"]
PLATFORMS = ['Reference', 'CPU', 'OpenCL', 'CUDA']
data_path = Path('tests')/'data'


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
        "Backbone": lambda: openawsem.functionTerms.basicTerms.con_term(oa),
        "Rama": lambda: openawsem.functionTerms.basicTerms.rama_term(oa),
        "Contact": lambda: openawsem.functionTerms.contactTerms.contact_term(oa),
        "Chain": lambda: openawsem.functionTerms.basicTerms.chain_term(oa),
        "Chi": lambda: openawsem.functionTerms.basicTerms.chi_term(oa),
        "Excluded": lambda: openawsem.functionTerms.basicTerms.excl_term(oa, periodic=False),
        "RamaProline": lambda: openawsem.functionTerms.basicTerms.rama_proline_term(oa),
        "RamaSSWeight": lambda: openawsem.functionTerms.basicTerms.rama_ssweight_term(oa, k_rama_ssweight=2*8.368, ssweight_file=data_path/f'{protein}-ssweight'),
        "Beta1": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_1(oa),
        "Beta2": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_2(oa),
        "Beta3": lambda: openawsem.functionTerms.hydrogenBondTerms.beta_term_3(oa),
        "Helical": lambda: openawsem.functionTerms.hydrogenBondTerms.helical_term(oa),
        "Pap1": lambda: openawsem.functionTerms.hydrogenBondTerms.pap_term_1(oa),
        "Pap2": lambda: openawsem.functionTerms.hydrogenBondTerms.pap_term_2(oa),
        "FragmentMemory": lambda: openawsem.functionTerms.templateTerms.fragment_memory_term(oa, frag_file_list_file=data_path/f'{protein}-single_frags.mem', npy_frag_table=data_path/f'{protein}-single_frags.npy', UseSavedFragTable=False),
        "DebyeHuckel": lambda: openawsem.functionTerms.debyeHuckelTerms.debye_huckel_term(oa, chargeFile=data_path/f'{protein}-charge.txt'),
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

    forceGroupTable = {"Backbone": 20, "Rama": 21, "Contact": 22, "Fragment": 23, "Membrane": 24, "ER": 25, "TBM_Q": 26, "Beta": 27, "Pap": 28, "Helical": 29,
                       "Q": 1, "Rg": 2, "Qc": 3, "Helix_orientation": 18, "Pulling": 19}

    termEnergies = pd.DataFrame(columns=["Step"] + COLUMNS)

    for step in range(len(pdb_trajectory)):
        simulation.context.setPositions(pdb_trajectory.openmm_positions(step))
        e = []
        for term in COLUMNS:
            g = {forceGroupTable[term]} if forceGroupTable[term] != -1 else -1
            state = simulation.context.getState(getEnergy=True, groups=g)
            termEnergy = state.getPotentialEnergy().value_in_unit(openawsem.unit_definitions.kilojoule_per_mole if term in COLUMNS else openawsem.unit_definitions.kilocalories_per_mole)
            e.append(termEnergy)
        termEnergies.loc[step] = [step] + e

    return termEnergies

def save_energies():
    simulation_platform = 'Reference'
    for protein in PROTEINS:
        termEnergies = analyze(protein, simulation_platform)
        termEnergies.to_csv(f'{protein}_energies.csv', index=False, float_format='%.6f')

def benchmark(protein, simulation_platform, timing_function = time_once, n_steps=100):
    chain = openawsem.helperFunctions.myFunctions.getAllChains(f"{protein}-crystal_structure.pdb")
    seq = openawsem.helperFunctions.myFunctions.read_fasta(f"{protein}-crystal_structure.fasta")
    pdb_trajectory = md.load(f'{protein}-movie.dcd', top=f"{protein}-openmmawsem.pdb")

    benchmark_data=[]
    for force_name in ['Backbone', 'Rama', 'Contact', 'Chain', 'Chi', 'Excluded', 'RamaProline', 'RamaSSWeight', 'Beta1', 'Beta2', 'Beta3', 'Helical', 'Pap1', 'Pap2', 'FragmentMemory', 'DebyeHuckel','All']:
        # Setup forces
        oa = openawsem.OpenMMAWSEMSystem(f"{protein}-openmmawsem.pdb",
                                         chains=chain,
                                         k_awsem=1.0,
                                         xml_filename=openawsem.xml,
                                         seqFromPdb=seq,
                                         includeLigands=False)
        
        forces, setup_time = timing_function(lambda: set_up_forces(oa, protein, force_name if force_name!='All' else None))
        oa.addForcesWithDefaultForceGroup(forces)
        print(f"{protein,simulation_platform,force_name} setup took {setup_time:.2f} ms")

        # Initialize simulation
        platform = openmm.Platform.getPlatformByName(simulation_platform)
        integrator = openmm.LangevinIntegrator(300*openawsem.unit_definitions.kelvin, 1/openawsem.unit_definitions.picosecond, .02*openawsem.unit_definitions.femtoseconds) # Timestep small to compute the force in the minimized step.
        simulation = openmm.app.Simulation(oa.pdb.topology, oa.system, integrator, platform)
        #simulation.context.setPositions(oa.pdb.positions)
        simulation.context.setPositions(pdb_trajectory.openmm_positions(0)) #Load positions after minimization

        # Run for n_steps and measure computation time
        _, simulation_time = timing_function(lambda: simulation.step(n_steps))
        print(f"{protein,simulation_platform,force_name} {n_steps} steps took {simulation_time:.2f} ms")

        benchmark_data.append([protein, simulation_platform, force_name, setup_time, simulation_time])
    return benchmark_data

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
        tolerance = 1e-5
        for protein in PROTEINS:
            calculated_energies = analyzed_data(protein, platform)
            saved_energies = pd.read_csv(data_path/f'{protein}_energies.csv')

            assert column in calculated_energies.columns, f"Column {column} not found in calculated energies for protein {protein} on platform {platform}"
            assert column in saved_energies.columns, f"Column {column} not found in saved energies for protein {protein} on platform {platform}"
            assert np.allclose(calculated_energies[column], saved_energies[column], atol=tolerance), f"Energy terms comparison failed for protein {protein} on column {column} on platform {platform}"


if __name__ == '__main__':
    # analyze("1brs", "Reference")
    data = []
    for protein in PROTEINS:
        for platform in PLATFORMS:
            data += benchmark(protein, platform, time_once)
    data = pd.DataFrame(data, columns=['protein', 'simulation_platform', 'force_name', 'setup_time', 'simulation_time'])
    data.to_csv('Benchmark_data_once.csv')

    data = []
    for protein in PROTEINS:
        for platform in PLATFORMS:
            data += benchmark(protein, platform, time_many)
    data = pd.DataFrame(data, columns=['protein', 'simulation_platform', 'force_name', 'setup_time', 'simulation_time'])
    data.to_csv('Benchmark_data_many.csv')