#!/usr/bin/env python3
"""
Regression test for membrane protein contact term with z_dependent=True.

This test reruns the 2xov trajectory and compares computed energies against
reference values from info.dat. The 2xov simulation was run with:
    contact_term(oa, k_contact=1*4.184, z_dependent=True, inMembrane=True, 
                 membrane_center=0, k_relative_mem=1, burialPartOn=True)

A bug in the energy string construction caused malformed expressions.
This test catches that by verifying energies match the reference.
"""
import pytest
import pandas as pd
import mdtraj as md
import numpy as np
import openmm
import openawsem
from pathlib import Path
from openmm.unit import angstrom, kilocalorie_per_mole

data_path = Path('tests') / 'data'


def load_reference_energies():
    """Load reference energies from info.dat."""
    # info.dat is whitespace-separated with header
    df = pd.read_csv(data_path / '2xov-info.dat', sep=r'\s+')
    return df


def setup_2xov_system_with_membrane_contact():
    """
    Set up OpenMMAWSEMSystem for 2xov with z_dependent contact term.
    
    This replicates the setup used in the original simulation (submode 3):
        contact_term(oa, k_contact=1*4.184, z_dependent=True, inMembrane=True,
                     membrane_center=0, k_relative_mem=1, burialPartOn=True)
    """
    chain = openawsem.helperFunctions.myFunctions.getAllChains(data_path / "2xov-crystal_structure.pdb")
    seq = openawsem.helperFunctions.myFunctions.read_fasta(data_path / "2xov-crystal_structure.fasta")
    
    oa = openawsem.OpenMMAWSEMSystem(
        data_path / "2xov-openmmawsem.pdb",
        chains=chain,
        k_awsem=1.0,
        xml_filename=openawsem.xml,
        seqFromPdb=seq,
        includeLigands=False
    )
    
    # Add contact term with z_dependent=True, inMembrane=True (submode 3)
    contact_force = openawsem.functionTerms.contactTerms.contact_term(
        oa,
        k_contact=1 * 4.184,  # 1 kcal/mol in kJ/mol
        z_dependent=True,
        inMembrane=True,
        membrane_center=0 * angstrom,
        k_relative_mem=1,
        parametersLocation=data_path,
        gammaName="2xov-gamma.dat",
        burialGammaName="2xov-burial_gamma.dat",
        membraneGammaName="2xov-membrane_gamma.dat",
        burialPartOn=True,
    )
    oa.system.addForce(contact_force)
    
    # Add membrane term
    membrane_force = openawsem.functionTerms.membraneTerms.membrane_preassigned_term(
        oa,
        k=5 * kilocalorie_per_mole,
        membrane_center=0 * angstrom,
        zimFile=str(data_path / "2xov-PredictedZim"),
    )
    oa.system.addForce(membrane_force)
    
    return oa


def create_simulation(oa, platform_name='Reference'):
    """Create a simulation from the OpenMMAWSEMSystem."""
    platform = openmm.Platform.getPlatformByName(platform_name)
    integrator = openmm.LangevinIntegrator(
        300 * openawsem.unit_definitions.kelvin,
        1 / openawsem.unit_definitions.picosecond,
        2 * openawsem.unit_definitions.femtoseconds
    )
    simulation = openmm.app.Simulation(oa.pdb.topology, oa.system, integrator, platform)
    return simulation


class TestMembraneProteinEnergies:
    """Test membrane protein energies by trajectory rerun against info.dat reference."""
    
    # Force group mapping (from test_energies.py)
    FORCE_GROUPS = {
        "Contact": 22,
        "Membrane": 24,
    }
    
    # info.dat stores energies in kcal/mol, but OpenMM uses kJ/mol
    KCAL_TO_KJ = 4.184
    
    def test_2xov_contact_energy_matches_reference(self):
        """
        Regression test: Contact energy from trajectory rerun must match info.dat.
        
        This test:
        1. Sets up 2xov with z_dependent=True, inMembrane=True contact term
        2. Reruns the saved trajectory (movie.dcd)
        3. Compares Contact energy against reference values from info.dat
        
        The bug in contact_term caused malformed energy strings, resulting in
        wrong energies or OpenMM parse errors.
        """
        # Load reference energies (in kcal/mol)
        ref = load_reference_energies()
        
        # Setup system
        oa = setup_2xov_system_with_membrane_contact()
        simulation = create_simulation(oa)
        
        # Load trajectory
        trajectory = md.load(data_path / '2xov-movie.dcd', top=data_path / "2xov-openmmawsem.pdb")
        
        # Tolerance for energy comparison (kJ/mol)
        tolerance = 1.0  # Allow 1 kJ/mol difference
        
        for step in range(len(trajectory)):
            simulation.context.setPositions(trajectory.openmm_positions(step))
            
            # Get Contact energy (force group 22) in kJ/mol
            state = simulation.context.getState(getEnergy=True, groups={self.FORCE_GROUPS["Contact"]})
            contact_energy = state.getPotentialEnergy().value_in_unit(openawsem.unit_definitions.kilojoule_per_mole)
            
            # Convert reference from kcal/mol to kJ/mol
            ref_contact = ref.loc[step, 'Contact'] * self.KCAL_TO_KJ
            
            assert np.isclose(contact_energy, ref_contact, atol=tolerance), \
                f"Step {step}: Contact energy mismatch. Got {contact_energy:.2f}, expected {ref_contact:.2f} kJ/mol"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
