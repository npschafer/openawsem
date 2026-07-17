import pytest
import shutil
from pathlib import Path
from openawsem.helperFunctions.convertOpenmmTrajectoryToStandardMovie import (
    get_seq_dic,
    convert_openMM_to_standard_pdb,
)

data_path = Path(__file__).parent / 'data' / 'DNA_reconstruction'


def parse_pdb_residues(pdb_text):
    """Parse PDB text and return a dict mapping (chain, resnum) -> resname."""
    residues = {}
    for line in pdb_text.splitlines():
        if line.startswith('ATOM') and len(line) > 26:
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26])
            residues[(chain, resnum)] = resname
    return residues


def assert_residues_converted(pdb_text):
    """Assert that all AWSEM placeholder residues were converted."""
    assert 'NGP' not in pdb_text
    assert 'IGL' not in pdb_text
    assert 'IPR' not in pdb_text


def assert_chain_boundary_residues(residues):
    """Assert correct residues at chain boundaries."""
    assert residues.get(('A', 1)) == 'MET'
    assert residues.get(('A', 274)) == 'ASP'
    assert residues.get(('B', 1)) == 'MET'
    assert residues.get(('B', 313)) == 'GLU'


class TestFixAminoacids:
    """Tests for the fix_aminoacids reconstruction workflow."""

    def test_reconstruction_with_multiple_chain_fasta(self, tmp_path):
        """Test reconstruction with properly formatted multi-chain FASTA."""
        pdb_copy = tmp_path / 'test.pdb'
        shutil.copy(data_path / 'bp40_frame5000_awsem.pdb', pdb_copy)
        
        seq_dic = get_seq_dic(fasta=data_path / 'protein_multiple.seq')
        convert_openMM_to_standard_pdb(fileName=pdb_copy, seq_dic=seq_dic, back=False)
        
        result = pdb_copy.read_text()
        assert_residues_converted(result)
        assert_chain_boundary_residues(parse_pdb_residues(result))

    def test_reconstruction_with_single_sequence_fasta(self, tmp_path):
        """Test reconstruction with deprecated single-sequence file (no chain headers)."""
        pdb_copy = tmp_path / 'test.pdb'
        shutil.copy(data_path / 'bp40_frame5000_awsem.pdb', pdb_copy)
        
        with pytest.warns(DeprecationWarning, match="No '>CRYSTAL_STRUCTURE:' chain headers found"):
            seq_dic = get_seq_dic(fasta=data_path / 'protein_single.seq')
        
        assert 'Unknown' in seq_dic
        assert len(seq_dic) == 1
        
        convert_openMM_to_standard_pdb(fileName=pdb_copy, seq_dic=seq_dic, back=False)
        
        result = pdb_copy.read_text()
        assert_residues_converted(result)
        assert_chain_boundary_residues(parse_pdb_residues(result))


class TestGetSeqDic:
    """Tests for the get_seq_dic function."""

    def test_multiple_chains(self):
        """Test parsing a valid multi-chain FASTA file."""
        seq_dic = get_seq_dic(fasta=data_path / 'protein_multiple.seq')

        assert len(seq_dic) == 2
        assert len(seq_dic['A']) == 274
        assert len(seq_dic['B']) == 313
        assert seq_dic['A'].startswith('MAYVEII')
        assert seq_dic['B'].startswith('MGPYLQI')

    def test_single_sequence_deprecated(self):
        """Test that plain sequence files emit deprecation warning."""
        with pytest.warns(DeprecationWarning, match="No '>CRYSTAL_STRUCTURE:' chain headers found"):
            seq_dic = get_seq_dic(fasta=data_path / 'protein_single.seq')
        
        assert 'Unknown' in seq_dic
        assert len(seq_dic) == 1
        assert len(seq_dic['Unknown']) == 587  # 274 + 313
        assert seq_dic['Unknown'].startswith('MAYVEII'), f"Sequence starts incorrectly: {seq_dic['Unknown'][:10]}"
        assert seq_dic['Unknown'].endswith('YYPE'), f"Sequence ends incorrectly: {seq_dic['Unknown'][-10:]}"
