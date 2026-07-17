"""
All-atom reconstruction from AWSEM+3SPN2 coarse-grained models.

Pipeline: CG PDB → fix_aminoacids (optional) → SCWRL4 (protein) + DNAbackmap (DNA) → merge → PDBFixer → OpenMM minimize

Requirements:
    - SCWRL4: http://dunbrack.fccc.edu/SCWRL4.php
    - DNAbackmap: GENESIS package (needs TCAG_fragment.txt in same dir)
    - OpenMM + PDBFixer: pip install openmm pdbfixer

Usage:
    awsem reconstruct input.pdb                         # Output: input_allatom.pdb
    awsem reconstruct input.pdb -f protein.seq          # Fix residue names first
    awsem reconstruct input.pdb -o output.pdb           # Custom output
    awsem reconstruct input.pdb --steps 10000           # More minimization
"""

import os
import sys
import subprocess
import shutil
import argparse
import tempfile
from pathlib import Path
from os import PathLike

# Tool paths - can be overridden via environment variables
SCWRL = os.environ.get('SCWRL4_PATH', '/usr/local/bin/scwrl4/Scwrl4')
DNABACKMAP = os.environ.get('DNABACKMAP_PATH', '/usr/local/bin/DNAbackmap')

# Standard amino acid residue names
PROTEIN_RES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

# AWSEM placeholder residue names (before fix_aminoacids converts them)
AWSEM_PROTEIN_RES = {'NGP', 'IGL', 'IPR'}

# All protein residue names (standard + AWSEM placeholders)
ALL_PROTEIN_RES = PROTEIN_RES | AWSEM_PROTEIN_RES

# DNA residue names (3SPN2 coarse-grained)
DNA_RES = {'DA', 'DC', 'DG', 'DT'}

# Backbone atoms to freeze during minimization
BACKBONE = {
    'N', 'CA', 'C', 'O', 'CB',  # Protein
    "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "P"  # DNA
}


def check_tool(path: str | PathLike, name: str) -> str | None:
    """Check if a tool exists and is executable.
    
    Args:
        path: Path to the tool executable
        name: Human-readable name for error messages
        
    Returns:
        Error message string if tool is missing/not executable, None if OK
    """
    tool_path = Path(path)
    if not tool_path.exists():
        return f"{name} not found at: {tool_path}"
    if not os.access(tool_path, os.X_OK):
        return f"{name} not executable: {tool_path}"
    return None


def check_tools(scwrl_path: str | None = None, 
                dnabackmap_path: str | None = None) -> dict[str, bool]:
    """Check which external tools are available.
    
    Args:
        scwrl_path: Override path to SCWRL4 executable
        dnabackmap_path: Override path to DNAbackmap executable
        
    Returns:
        Dict mapping tool name to availability (True = available)
    """
    scwrl = scwrl_path or SCWRL
    dnabackmap = dnabackmap_path or DNABACKMAP
    
    status = {
        'scwrl': check_tool(scwrl, 'SCWRL4') is None,
        'dnabackmap': check_tool(dnabackmap, 'DNAbackmap') is None
    }
    
    # DNAbackmap also needs its fragment library
    if status['dnabackmap']:
        lib = Path(dnabackmap).parent / 'TCAG_fragment.txt'
        if not lib.exists():
            status['dnabackmap'] = False
    
    return status


def fix_aminoacids(pdb_in: str | PathLike, seq_file: str | PathLike, 
                   pdb_out: str | PathLike) -> bool:
    """Fix residue names using sequence file.
    
    Uses internal AWSEM functions to convert placeholder residue names
    (NGP, IGL, IPR) to standard amino acid names based on the sequence.
    
    Args:
        pdb_in: Input CG PDB file
        seq_file: Sequence file (FASTA format)
        pdb_out: Output PDB file with fixed residue names
        
    Returns:
        True if successful, False otherwise
    """
    try:
        from openawsem.helperFunctions.convertOpenmmTrajectoryToStandardMovie import (
            get_seq_dic, convert_openMM_to_standard_pdb
        )
        
        # Copy input to output first (convert_openMM_to_standard_pdb modifies in place)
        shutil.copy(pdb_in, pdb_out)
        
        # Get sequence dictionary from FASTA file
        seq_dic = get_seq_dic(fasta=seq_file)
        
        # Convert residue names in place
        convert_openMM_to_standard_pdb(fileName=pdb_out, seq_dic=seq_dic, back=False)
        
        pdb_out = Path(pdb_out)
        return pdb_out.exists() and pdb_out.stat().st_size > 0
        
    except Exception as e:
        print(f"Warning: fix_aminoacids error: {e}")
        return False


def extract_protein(pdb_in: str | PathLike, pdb_out: str | PathLike) -> int:
    """Extract protein atoms (AWSEM backbone: CA,C,O,N,H,CB) from CG PDB.
    
    Recognizes both standard amino acid names and AWSEM placeholder names (NGP, IGL, IPR).
    
    Args:
        pdb_in: Input CG PDB file
        pdb_out: Output PDB file containing only protein atoms
        
    Returns:
        Number of protein atoms extracted
    """
    lines = []
    with open(pdb_in) as f:
        for line in f:
            if line[:4] == 'ATOM' and line[17:20].strip() in ALL_PROTEIN_RES:
                lines.append(line)
    
    with open(pdb_out, 'w') as f:
        f.write(''.join(lines) + 'END\n')
    
    return len(lines)


def extract_dna(pdb_in: str | PathLike, pdb_out: str | PathLike) -> int:
    """Extract DNA atoms (3SPN2: Base,Sugar,Phosphate) and rename for DNAbackmap.
    
    Renames atom types: A/T/G/C -> DB, S -> DS, P -> DP
    
    Args:
        pdb_in: Input CG PDB file
        pdb_out: Output PDB file with renamed DNA atoms
        
    Returns:
        Number of DNA atoms extracted
    """
    # Atom name mapping: 3SPN2 -> DNAbackmap
    atom_map = {
        'A': 'DB  ', 'T': 'DB  ', 'G': 'DB  ', 'C': 'DB  ',  # Base
        'S': 'DS  ',  # Sugar
        'P': 'DP  '   # Phosphate
    }
    
    lines = []
    with open(pdb_in) as f:
        for line in f:
            if line[:4] == 'ATOM' and line[17:20].strip() in DNA_RES:
                atom_name = line[12:16].strip()
                new_name = atom_map.get(atom_name)
                if new_name:
                    lines.append(line[:12] + new_name + line[16:])
    
    with open(pdb_out, 'w') as f:
        f.write(''.join(lines) + 'END\n')
    
    return len(lines)


def run_scwrl(pdb_in: str | PathLike, pdb_out: str | PathLike,
              scwrl_path: str | None = None) -> bool:
    """Add side chains to protein backbone using SCWRL4.
    
    Args:
        pdb_in: Input PDB file with protein backbone
        pdb_out: Output PDB file with added side chains
        scwrl_path: Override path to SCWRL4 executable
        
    Returns:
        True if successful, False otherwise
    """
    scwrl = scwrl_path or SCWRL
    try:
        result = subprocess.run(
            [scwrl, '-i', str(pdb_in), '-o', str(pdb_out)],
            capture_output=True, text=True, timeout=600
        )
        if result.returncode != 0:
            print(f"SCWRL4 error: {result.stderr.strip() or result.stdout.strip()}")
            return False
        out_path = Path(pdb_out)
        if not (out_path.exists() and out_path.stat().st_size > 0):
            print("SCWRL4 error: output file not created")
            return False
        return True
    except subprocess.TimeoutExpired:
        print("SCWRL4 error: timed out after 600s")
        return False
    except Exception as e:
        print(f"SCWRL4 error: {e}")
        return False


def run_dnabackmap(pdb_in: str | PathLike, pdb_out: str | PathLike,
                   dnabackmap_path: str | None = None,
                   work_dir: str | PathLike | None = None) -> bool:
    """Reconstruct all-atom DNA from CG using DNAbackmap.
    
    Note: DNAbackmap must run from its install directory where TCAG_fragment.txt is located.
    
    Args:
        pdb_in: Input CG DNA PDB file (with DB/DS/DP atom names)
        pdb_out: Output all-atom DNA PDB file
        dnabackmap_path: Override path to DNAbackmap executable
        work_dir: Working directory for temp files (avoids race conditions)
        
    Returns:
        True if successful, False otherwise
    """
    dnabackmap = Path(dnabackmap_path or DNABACKMAP)
    dnabackmap_dir = dnabackmap.parent
    
    # Use provided work_dir or create unique temp files in tool directory
    if work_dir:
        work = Path(work_dir)
    else:
        # Fallback: use tool directory with PID to avoid collisions
        work = dnabackmap_dir
    
    # Use PID for unique temp file names to avoid race conditions
    pid = os.getpid()
    cg_file = f'temp_cg_{pid}.pdb'
    aa_file = f'temp_aa_{pid}.pdb'
    input_file = f'input_temp_{pid}.txt'
    
    try:
        # Copy input to DNAbackmap directory (required by tool)
        shutil.copy(pdb_in, dnabackmap_dir / cg_file)
        
        # Create input file
        (dnabackmap_dir / input_file).write_text(f'FILENAME {cg_file} {aa_file}\n')
        
        result = subprocess.run(
            [str(dnabackmap), input_file],
            cwd=dnabackmap_dir,
            capture_output=True, text=True, timeout=300
        )
        
        if result.returncode != 0:
            print(f"DNAbackmap error: {result.stderr.strip() or result.stdout.strip()}")
            return False
        
        aa_path = dnabackmap_dir / aa_file
        if aa_path.exists():
            shutil.copy(aa_path, pdb_out)
            return True
        
        print("DNAbackmap error: output file not created")
        return False
        
    except subprocess.TimeoutExpired:
        print("DNAbackmap error: timed out after 300s")
        return False
    except Exception as e:
        print(f"DNAbackmap error: {e}")
        return False
    finally:
        # Clean up temp files
        for f in [cg_file, aa_file, input_file]:
            try:
                (dnabackmap_dir / f).unlink()
            except OSError:
                pass


def merge(prot_pdb: str | PathLike, dna_pdb: str | PathLike, 
          out_pdb: str | PathLike) -> None:
    """Merge protein and DNA PDBs with renumbered atoms.
    
    Args:
        prot_pdb: Input protein PDB file
        dna_pdb: Input DNA PDB file
        out_pdb: Output merged PDB file
    """
    lines = []
    atom_num = 1
    
    for pdb_file in [prot_pdb, dna_pdb]:
        with open(pdb_file) as f:
            for line in f:
                if line[:4] == 'ATOM':
                    # Renumber atom
                    lines.append(f'{line[:6]}{atom_num:5d}{line[11:]}')
                    atom_num += 1
    
    with open(out_pdb, 'w') as f:
        f.write(''.join(lines) + 'END\n')


def fix_minimize(pdb_in: str | PathLike, pdb_out: str | PathLike, 
                 steps: int = 5000, freeze_backbone: bool = True) -> bool:
    """Fix missing atoms/hydrogens with PDBFixer, then energy minimize with OpenMM.
    
    Args:
        pdb_in: Input PDB file
        pdb_out: Output minimized PDB file
        steps: Maximum minimization iterations
        freeze_backbone: Whether to freeze backbone atoms during minimization
        
    Returns:
        True if successful, False otherwise
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile, Modeller, ForceField, NoCutoff, Simulation
    from openmm import LangevinMiddleIntegrator, Platform, unit
    
    # Fix structure with PDBFixer
    fixer = PDBFixer(filename=str(pdb_in))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    # Set up OpenMM system
    ff = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    modeller = Modeller(fixer.topology, fixer.positions)
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff, removeCMMotion=True)
    
    # Freeze backbone atoms by setting mass to 0
    if freeze_backbone:
        for atom in modeller.topology.atoms():
            if atom.name in BACKBONE:
                system.setParticleMass(atom.index, 0 * unit.dalton)
    
    # Set up integrator and simulation
    integrator = LangevinMiddleIntegrator(50 * unit.kelvin, 1 / unit.picosecond, 1 * unit.femtoseconds)
    
    try:
        platform = Platform.getPlatformByName("CUDA")
    except Exception:
        platform = Platform.getPlatformByName("CPU")
    
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize energy
    simulation.minimizeEnergy(maxIterations=steps)
    
    # Write output
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(pdb_out, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    
    print(f"Energy: {state.getPotentialEnergy()}")
    return True


def reconstruct(input_pdb: str | PathLike, output_pdb: str | PathLike | None = None,
                seq_file: str | PathLike | None = None, steps: int = 5000,
                freeze_backbone: bool = True,
                scwrl_path: str | None = None,
                dnabackmap_path: str | None = None) -> Path:
    """Reconstruct all-atom structure from AWSEM+3SPN2 coarse-grained model.
    
    Args:
        input_pdb: Input CG PDB file
        output_pdb: Output all-atom PDB file (default: input_allatom.pdb)
        seq_file: Sequence file for fix_aminoacids (optional)
        steps: Minimization steps
        freeze_backbone: Whether to freeze backbone during minimization
        scwrl_path: Path to SCWRL4 executable (optional)
        dnabackmap_path: Path to DNAbackmap executable (optional)
        
    Returns:
        Path to output PDB file
    """
    inp = Path(input_pdb)
    
    if output_pdb is None:
        stem = inp.stem.replace('_awsem', '').replace('_fixed', '')
        out = inp.parent / f"{stem}_allatom.pdb"
    else:
        out = Path(output_pdb)
    
    # Create temp directory
    tmp = Path(tempfile.mkdtemp())
    
    try:
        current_pdb = inp
        
        # Fix residue names if sequence file provided
        if seq_file:
            fixed_pdb = tmp / 'fixed.pdb'
            print(f"Running fix_aminoacids with {seq_file}...")
            if fix_aminoacids(current_pdb, seq_file, fixed_pdb):
                current_pdb = fixed_pdb
                print("Residue names fixed successfully")
            else:
                print("Warning: Continuing without fix_aminoacids")
        
        # Extract protein and DNA
        prot_cg = tmp / 'prot_cg.pdb'
        prot_aa = tmp / 'prot_aa.pdb'
        dna_cg = tmp / 'dna_cg.pdb'
        dna_aa = tmp / 'dna_aa.pdb'
        merged = tmp / 'merged.pdb'
        
        print(f"Input: {inp}\nOutput: {out}")
        
        n_prot = extract_protein(current_pdb, prot_cg)
        n_dna = extract_dna(current_pdb, dna_cg)
        print(f"Extracted: {n_prot} protein atoms, {n_dna} DNA atoms")
        
        # Check which tools are available
        tools = check_tools(scwrl_path=scwrl_path, dnabackmap_path=dnabackmap_path)
        
        # Reconstruct protein sidechains with SCWRL4
        prot_reconstructed = False
        if n_prot > 0:
            if tools['scwrl']:
                print("Running SCWRL4...")
                if run_scwrl(prot_cg, prot_aa, scwrl_path=scwrl_path):
                    prot_reconstructed = True
                else:
                    print("Warning: SCWRL4 failed, using backbone-only protein")
            else:
                print("Warning: SCWRL4 not available, skipping side-chain reconstruction")
        
        # Reconstruct all-atom DNA with DNAbackmap
        dna_reconstructed = False
        if n_dna > 0:
            if tools['dnabackmap']:
                print("Running DNAbackmap...")
                if run_dnabackmap(dna_cg, dna_aa, dnabackmap_path=dnabackmap_path, work_dir=tmp):
                    dna_reconstructed = True
                else:
                    print("Warning: DNAbackmap failed, using CG DNA")
            else:
                print("Warning: DNAbackmap not available, skipping DNA reconstruction")
        
        # Determine what files to use for merging
        prot_file = prot_aa if prot_reconstructed else (prot_cg if n_prot > 0 else None)
        dna_file = dna_aa if dna_reconstructed else (dna_cg if n_dna > 0 else None)
        
        # Merge reconstructed structures first (preserves chain IDs)
        print("Merging structures...")
        if prot_file and dna_file:
            merge(prot_file, dna_file, merged)
        elif prot_file:
            shutil.copy(prot_file, merged)
        elif dna_file:
            shutil.copy(dna_file, merged)
        else:
            raise RuntimeError("No protein or DNA atoms found in input")
        
        # Run PDBFixer and minimize on merged structure
        # (Only works if both protein and DNA are all-atom)
        if prot_reconstructed or dna_reconstructed:
            can_minimize = True
            # If DNA wasn't reconstructed, CG residues will fail in force field
            if n_dna > 0 and not dna_reconstructed:
                print("Warning: Skipping minimization (CG DNA not supported by force field)")
                can_minimize = False
            # If protein wasn't reconstructed, CG residues may fail
            if n_prot > 0 and not prot_reconstructed:
                print("Warning: Skipping minimization (CG protein not fully supported)")
                can_minimize = False
            
            if can_minimize:
                print("Running PDBFixer and minimizing merged structure...")
                fix_minimize(merged, out, steps, freeze_backbone)
            else:
                shutil.copy(merged, out)
        else:
            print("Warning: No reconstruction tools available, output contains CG structure")
            shutil.copy(merged, out)
        
        print(f"Done: {out}")
        
        return out
        
    finally:
        # Clean up temp directory
        shutil.rmtree(tmp, ignore_errors=True)


def main(args=None):
    """CLI entry point for reconstruction."""
    parser = argparse.ArgumentParser(
        description='Reconstruct all-atom structure from AWSEM+3SPN2 coarse-grained model'
    )
    parser.add_argument('input', help='Input CG PDB file')
    parser.add_argument('-o', '--output', help='Output all-atom PDB file')
    parser.add_argument(
        '-f', '--seq', 
        help='Sequence file (FASTA) for fix_aminoacids (e.g., protein.seq)'
    )
    parser.add_argument(
        '--steps', type=int, default=5000,
        help='Minimization steps (default: 5000)'
    )
    parser.add_argument(
        '--no-freeze-backbone', action='store_true',
        help='Do not freeze backbone atoms during minimization'
    )
    parser.add_argument(
        '--scwrl', metavar='PATH',
        help='Path to SCWRL4 executable (default: $SCWRL4_PATH or /usr/local/bin/scwrl4/Scwrl4)'
    )
    parser.add_argument(
        '--dnabackmap', metavar='PATH',
        help='Path to DNAbackmap executable (default: $DNABACKMAP_PATH or /usr/local/bin/DNAbackmap)'
    )
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    try:
        reconstruct(
            input_pdb=args.input,
            output_pdb=args.output,
            seq_file=args.seq,
            steps=args.steps,
            freeze_backbone=not args.no_freeze_backbone,
            scwrl_path=args.scwrl,
            dnabackmap_path=args.dnabackmap
        )
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
