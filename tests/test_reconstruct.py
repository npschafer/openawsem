"""Tests for the atomistic reconstruction workflow.

These tests focus on:
1. Core logic that CAN break (extraction, renaming, merging)
2. Integration with real test data
3. Error handling that users will actually encounter

Not tested (by design):
- Trivial wrappers (check_tool is just os.path.exists + os.access)
- Heavily mocked subprocess calls that don't test real behavior
- Constants (if PROTEIN_RES is wrong, fix the code, not the test)
"""
import pytest
import subprocess
from pathlib import Path
from unittest.mock import patch

# Test data path
data_path = Path(__file__).parent / 'data' / 'DNA_reconstruction'


class TestExtractProtein:
    """Tests protein extraction - the core logic that separates protein from DNA."""

    def test_extracts_only_protein_from_mixed_pdb(self, tmp_path):
        """Extract protein atoms from CG PDB, excluding all DNA."""
        from openawsem.helperFunctions.reconstruct import extract_protein, ALL_PROTEIN_RES
        
        pdb_in = data_path / 'bp40_frame5000_awsem.pdb'
        pdb_out = tmp_path / 'protein.pdb'
        
        count = extract_protein(pdb_in, pdb_out)
        
        assert count > 0, "Test data should contain protein"
        
        # Verify ONLY protein residues (this catches bugs in the residue filter)
        content = pdb_out.read_text()
        for line in content.splitlines():
            if line.startswith('ATOM'):
                resname = line[17:20].strip()
                assert resname in ALL_PROTEIN_RES, f"Leaked non-protein residue: {resname}"
        
        # DNA should never appear (catches bugs where filter is too permissive)
        for dna_res in ['DA ', 'DC ', 'DG ', 'DT ']:
            assert dna_res not in content, f"DNA residue {dna_res.strip()} leaked into protein output"


class TestExtractDNA:
    """Tests DNA extraction and critical atom renaming for DNAbackmap."""

    def test_extracts_and_renames_dna_atoms(self, tmp_path):
        """Extract DNA and verify atom renaming (A/T/G/C->DB, S->DS, P->DP).
        
        This is critical because DNAbackmap expects specific atom names.
        Wrong names = silent failure in reconstruction.
        """
        from openawsem.helperFunctions.reconstruct import extract_dna
        
        # Create synthetic 3SPN2 DNA to test all atom types
        cg_pdb = tmp_path / 'cg_dna.pdb'
        cg_pdb.write_text(
            "ATOM      1  A   DA  C   1       0.000   0.000   0.000  1.00  0.00\n"
            "ATOM      2  S   DA  C   1       1.000   0.000   0.000  1.00  0.00\n"
            "ATOM      3  P   DA  C   2       2.000   0.000   0.000  1.00  0.00\n"
            "ATOM      4  T   DT  C   2       3.000   0.000   0.000  1.00  0.00\n"
            "ATOM      5  G   DG  C   3       4.000   0.000   0.000  1.00  0.00\n"
            "ATOM      6  C   DC  C   4       5.000   0.000   0.000  1.00  0.00\n"
            "END\n"
        )
        pdb_out = tmp_path / 'dna.pdb'
        
        count = extract_dna(cg_pdb, pdb_out)
        
        assert count == 6, "Should extract all 6 atoms"
        
        lines = [l for l in pdb_out.read_text().splitlines() if l.startswith('ATOM')]
        atom_names = [l[12:16].strip() for l in lines]
        
        # Exact expected output - catches any renaming bugs
        assert atom_names == ['DB', 'DS', 'DP', 'DB', 'DB', 'DB'], \
            f"Atom renaming failed: got {atom_names}"


class TestMerge:
    """Tests correct merging and atom renumbering."""

    def test_preserves_chain_ids_and_renumbers_atoms(self, tmp_path):
        """Merge must preserve chain IDs but renumber atoms sequentially.
        
        Chain ID preservation is critical for downstream analysis.
        """
        from openawsem.helperFunctions.reconstruct import merge
        
        prot_pdb = tmp_path / 'protein.pdb'
        dna_pdb = tmp_path / 'dna.pdb'
        merged_pdb = tmp_path / 'merged.pdb'
        
        # Deliberately use different chain IDs and non-sequential atom numbers
        prot_pdb.write_text(
            "ATOM    100  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
            "ATOM    101  C   ALA B   1       1.000   0.000   0.000  1.00  0.00\n"
            "END\n"
        )
        dna_pdb.write_text(
            "ATOM    500  P   DA  C   1       5.000   0.000   0.000  1.00  0.00\n"
            "ATOM    501  C5' DA  D   1       6.000   0.000   0.000  1.00  0.00\n"
            "END\n"
        )
        
        merge(prot_pdb, dna_pdb, merged_pdb)
        
        lines = [l for l in merged_pdb.read_text().splitlines() if l.startswith('ATOM')]
        
        # Atoms renumbered 1-4
        atom_nums = [int(l[6:11].strip()) for l in lines]
        assert atom_nums == [1, 2, 3, 4], f"Atom renumbering failed: {atom_nums}"
        
        # Chain IDs preserved (A, B, C, D)
        chain_ids = [l[21] for l in lines]
        assert chain_ids == ['A', 'B', 'C', 'D'], f"Chain IDs changed: {chain_ids}"


class TestFixAminoacids:
    """Tests residue name conversion using AWSEM internal functions."""

    def test_converts_placeholder_residues(self, tmp_path):
        """Convert NGP/IGL/IPR placeholders to standard amino acids.
        
        This is required before SCWRL4 can add sidechains.
        """
        from openawsem.helperFunctions.reconstruct import fix_aminoacids
        
        pdb_in = data_path / 'bp40_frame5000_awsem.pdb'
        pdb_out = tmp_path / 'fixed.pdb'
        seq_file = data_path / 'protein_multiple.seq'
        
        # First verify input actually has placeholders
        input_content = pdb_in.read_text()
        has_placeholders = any(p in input_content for p in ['NGP', 'IGL', 'IPR'])
        assert has_placeholders, "Test data should contain placeholder residues"
        
        result = fix_aminoacids(pdb_in, seq_file, pdb_out)
        
        assert result is True
        
        # Verify placeholders are gone
        output_content = pdb_out.read_text()
        for placeholder in ['NGP', 'IGL', 'IPR']:
            assert placeholder not in output_content, f"{placeholder} not converted"

    def test_returns_false_on_missing_sequence(self, tmp_path):
        """Gracefully handle missing sequence file instead of crashing."""
        from openawsem.helperFunctions.reconstruct import fix_aminoacids
        
        pdb_in = data_path / 'bp40_frame5000_awsem.pdb'
        pdb_out = tmp_path / 'fixed.pdb'
        
        result = fix_aminoacids(pdb_in, tmp_path / 'nonexistent.seq', pdb_out)
        
        assert result is False


class TestFixMinimize:
    """Tests PDBFixer + OpenMM minimization on valid structures."""

    def test_minimizes_valid_peptide(self, tmp_path):
        """Run minimization on a valid peptide structure.
        
        This catches force field compatibility issues and broken minimization logic.
        """
        pytest.importorskip('openmm', reason='OpenMM required')
        from openawsem.helperFunctions.reconstruct import fix_minimize
        
        pdb_in = tmp_path / 'input.pdb'
        pdb_out = tmp_path / 'output.pdb'
        
        # 3-residue ALA peptide with reasonable geometry
        pdb_in.write_text("""\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.382   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.768   1.215  1.00  0.00           C
ATOM      6  N   ALA A   2       3.315   1.550   0.000  1.00  0.00           N
ATOM      7  CA  ALA A   2       3.950   2.850   0.000  1.00  0.00           C
ATOM      8  C   ALA A   2       5.460   2.750   0.000  1.00  0.00           C
ATOM      9  O   ALA A   2       6.100   1.700   0.000  1.00  0.00           O
ATOM     10  CB  ALA A   2       3.500   3.700   1.200  1.00  0.00           C
ATOM     11  N   ALA A   3       6.000   3.950   0.000  1.00  0.00           N
ATOM     12  CA  ALA A   3       7.440   4.150   0.000  1.00  0.00           C
ATOM     13  C   ALA A   3       8.100   2.780   0.000  1.00  0.00           C
ATOM     14  O   ALA A   3       7.400   1.780   0.000  1.00  0.00           O
ATOM     15  CB  ALA A   3       7.900   4.950   1.200  1.00  0.00           C
END
""")
        
        fix_minimize(pdb_in, pdb_out, steps=10, freeze_backbone=False)
        
        assert pdb_out.exists()
        # Output should have more atoms (hydrogens added)
        n_out = sum(1 for l in pdb_out.read_text().splitlines() if l.startswith('ATOM'))
        assert n_out > 15, "PDBFixer should add hydrogens"


class TestIntegration:
    """End-to-end tests using real test data."""

    def test_full_extraction_merge_pipeline(self, tmp_path):
        """Extract protein+DNA from test data, merge, verify structure.
        
        This is the core pipeline that runs before external tools.
        """
        from openawsem.helperFunctions.reconstruct import extract_protein, extract_dna, merge
        
        pdb_in = data_path / 'bp40_frame5000_awsem.pdb'
        prot_pdb = tmp_path / 'protein.pdb'
        dna_pdb = tmp_path / 'dna.pdb'
        merged_pdb = tmp_path / 'merged.pdb'
        
        n_prot = extract_protein(pdb_in, prot_pdb)
        n_dna = extract_dna(pdb_in, dna_pdb)
        
        assert n_prot > 100, f"Expected many protein atoms, got {n_prot}"
        assert n_dna > 100, f"Expected many DNA atoms, got {n_dna}"
        
        merge(prot_pdb, dna_pdb, merged_pdb)
        
        # Verify merged atom count
        n_merged = sum(1 for l in merged_pdb.read_text().splitlines() if l.startswith('ATOM'))
        assert n_merged == n_prot + n_dna, "Merge lost atoms"
        
        # Verify we have expected chain IDs (A, B for protein, C, D for DNA)
        chains = set()
        for line in merged_pdb.read_text().splitlines():
            if line.startswith('ATOM'):
                chains.add(line[21])
        assert len(chains) == 4, f"Expected 4 chains (A,B,C,D), got {chains}"

    def test_cli_runs(self):
        """Verify CLI can be invoked (smoke test)."""
        result = subprocess.run(
            ['python', '-m', 'openawsem.helperFunctions.reconstruct', '--help'],
            capture_output=True, text=True
        )
        assert result.returncode == 0

    def test_protein_only_does_not_call_dnabackmap(self, tmp_path):
        """Protein-only input should never invoke DNAbackmap.
        
        This is critical: users without DNAbackmap installed should still
        be able to reconstruct protein-only AWSEM structures.
        """
        from openawsem.helperFunctions.reconstruct import (
            extract_protein, extract_dna, reconstruct
        )
        
        # Create protein-only PDB from test data
        pdb_in = data_path / 'bp40_frame5000_awsem.pdb'
        protein_only = tmp_path / 'protein_only.pdb'
        n_prot = extract_protein(pdb_in, protein_only)
        assert n_prot > 0
        
        # Verify no DNA in our test file
        n_dna = extract_dna(protein_only, tmp_path / 'dna_check.pdb')
        assert n_dna == 0, "Test file should have no DNA"
        
        # Mock run_dnabackmap to track if it's called
        with patch('openawsem.helperFunctions.reconstruct.run_dnabackmap') as mock_dnabackmap:
            with patch('openawsem.helperFunctions.reconstruct.run_scwrl') as mock_scwrl:
                # Make SCWRL4 "succeed" by copying input to output
                def fake_scwrl(pdb_in, pdb_out, scwrl_path=None):
                    import shutil
                    shutil.copy(pdb_in, pdb_out)
                    return True
                mock_scwrl.side_effect = fake_scwrl
                
                with patch('openawsem.helperFunctions.reconstruct.fix_minimize') as mock_minimize:
                    def fake_minimize(pdb_in, pdb_out, steps=5000, freeze_backbone=True):
                        import shutil
                        shutil.copy(pdb_in, pdb_out)
                        return True
                    mock_minimize.side_effect = fake_minimize
                    
                    out = reconstruct(protein_only, tmp_path / 'output.pdb')
        
        # DNAbackmap should NEVER be called for protein-only input
        mock_dnabackmap.assert_not_called()
        
        # SCWRL4 should be called (for protein)
        mock_scwrl.assert_called_once()
        
        # Output should exist
        assert out.exists()
