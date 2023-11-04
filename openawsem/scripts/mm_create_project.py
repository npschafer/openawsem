#!/usr/bin/env python3
import os
import argparse
import sys
import openawsem
import pandas as pd
from pathlib import Path
import contextlib
import subprocess


__location__ = openawsem.__location__
__author__ = 'Wei Lu'

def parse_arguments():
    """
    This function defines and parses command-line arguments for a protein simulation project template creation script.
    It uses the argparse module to define the expected arguments and their corresponding help messages.
    """
    
    # Create an argument parser with a description for the script
    parser = argparse.ArgumentParser(
        description=f"This Python 3 script automatically creates a protein simulation project template quickly and efficiently. Written by {__author__}"
    )

    # Define the expected arguments and their help messages
    parser.add_argument("proteins", nargs="*", help="Provide the names of the proteins (e.g., 1r69) or the target PDB files for the simulation, separated by spaces.")
    parser.add_argument("-c", "--chain", default="-1", help="Specify the chains to be simulated (e.g., 'ABC').")
    parser.add_argument("-d", "--debug", action="store_true", default=False, help="Enable debug mode.")
    parser.add_argument("--frag", action="store_true", default=False, help="Generate fragment memories.")
    parser.add_argument("--extended", action="store_true", default=False, help="Start from an extended structure generated using PyMOL (ensure it's installed). Supports single chain only.")
    parser.add_argument("--membrane", action="store_true", default=False, help="Enable membrane protein simulations.")
    parser.add_argument("--hybrid", action="store_true", default=False, help="Enable hybrid simulations.")
    parser.add_argument("--verbose", action="store_true", default=False, help="Enable verbose output.")
    parser.add_argument("--predict_ssweight_from_fasta", action="store_true", default=False, help="Predict secondary structure weight from FASTA sequence.")
    parser.add_argument("--keepIds", action="store_true", default=False, help="Preserve chain and residue index. By default, chains will be renamed from 'A' and indices will start from 1.")
    parser.add_argument("--keepLigands", action="store_true", default=False, help="Preserve ligands in the protein structure.")
    parser.add_argument("--test", action="store_true", default=False, help="Tests the current module")

      # Create a subparser for frag-related arguments
    frag_parser = parser.add_argument_group("frag", "Arguments for fragment memory generation. Only used if --frag is specified")
    frag_parser.add_argument("--frag_database", default=openawsem.data_path.blast, help="Specify the database for fragment generation.")
    frag_parser.add_argument("--frag_fasta", default=None, help="Provide the FASTA file for fragment generation.")
    frag_parser.add_argument("--frag_N_mem", type=int, default=20, help="Number of memories to generate per fragment.")
    frag_parser.add_argument("--frag_brain_damage", type=float, choices=[0, 0.5, 1, 2], default=0, help="Control the inclusion or exclusion of homologous protein structures for generating fragment memories.\n 0: include all homologs, 0.5: include only self-structures, 1: exclude all homologs, 2: include only non-homologous structures.")
    frag_parser.add_argument("--frag_fragmentLength", type=int, default=9, help="SLength of the fragments to be generated.")


    # Parse and return the command-line arguments
    return parser.parse_args()

class AWSEMSimulationProject:
    def __init__(self,args):
        self.data_path=__location__
        self.base_folder=Path.cwd() #Project folder
        self.args=args
                
    def run_command(self, command, stdout=None):

        if self.args.debug:
            print(' '.join(command))
        else:
            subprocess.run(command, check=True, shell=False, stdout=(open(stdout, "w") if stdout else None))

    @contextlib.contextmanager
    def change_directory(self,path):
        old_path = os.getcwd()
        try:
            os.chdir(path)
            yield
        finally:
            os.chdir(old_path)

            
    def log_commandline_args(self, log_file = 'create_project_commandline_args.txt'):
        with open(log_file, 'w') as f:
            f.write(' '.join(sys.argv))
            f.write('\n')

    def prepare_input_files_from_pdb(self, parent_folder="."):
        """
        Prepare input files from a PDB file.
        
        This function checks if the PDB file exists, extracts the name of the protein,
        and copies the PDB file to the 'original_pdbs' directory.
        
        Returns:
            name (str): Protein name without the file extension.
            pdb (str): PDB filename.
        """
        protein_path = Path(self.args.protein)
        parent_folder = Path(parent_folder)
        assert protein_path.exists(), f"The protein path {str(protein_path)} does not exist"
        
        #Get the name of the pdb
        name = protein_path.stem
        pdb = protein_path.name
        
        #Create a pdb backup
        original_pdbs_dir = parent_folder / "original_pdbs"
        original_pdbs_dir.mkdir(parents=True, exist_ok=True)
        (original_pdbs_dir / pdb).write_bytes(protein_path.read_bytes())
        
        return name, pdb

    def prepare_input_files_from_fasta(self, parent_folder="."):
        """
        Prepare input files from a FASTA file.
        
        This function creates a simulation folder from a FASTA file, extracts the name of the protein,
        generates a PDB file, and copies both the FASTA and PDB files to their respective directories.
        
        Returns:
            name (str): Protein name without the file extension.
            pdb (str): PDB filename.
        """
        print("Creating simulation folder from FASTA file.")
        protein_path = Path(self.args.protein)
        name = protein_path.stem

        try:
            self.run_command(["python3", __location__/"helperFunctions"/"fasta2pdb.py", name, "-f", self.args.protein])
            openawsem.helperFunctions.add_chain_to_pymol_pdb(f"{name}.pdb")
        except Exception as e:
            print(f"ERROR: Failed to convert FASTA to PDB. Exception: {e}")
            exit()

        original_fasta_dir = Path(parent_folder) / "original_fasta"
        original_pdbs_dir = Path(parent_folder) / "original_pdbs"
        original_fasta_dir.mkdir(parents=True, exist_ok=True)
        original_pdbs_dir.mkdir(parents=True, exist_ok=True)
        fasta_file_path = (original_fasta_dir / protein_path.name)
        fasta_file_path.write_bytes(protein_path.read_bytes())
        pdb_path = Path(f"{name}.pdb")
        crystal_structure_path=(Path(parent_folder) / "crystal_structure.pdb")
        crystal_structure_path.write_bytes(pdb_path.read_bytes())
        pdb = pdb_path.name
        
        return name, pdb

    def prepare_input_files_from_name(self, parent_folder="."):
        """
        Prepare input files from a protein name.
        
        This function downloads the PDB file for the given protein name.
        
        Returns:
            name (str): Protein name without the file extension.
            pdb (str): PDB filename.
        """
        name = self.args.protein
        pdb = f"{name}.pdb"
        pdb_list = [name]
        parent_folder = Path(parent_folder)
        
        try:
            openawsem.helperFunctions.downloadPdb(pdb_list, location=parent_folder/'original_pdbs')
        except Exception as e:
            print(f"ERROR: Failed to download PDB file. Exception: {e}")
            exit()

        return name, pdb

    def process_pdb_files(self):
        """
        Process the PDB files by cleaning, preparing, and generating additional required files.
        """
        removeHeterogens = False if self.args.keepLigands is True else True
        chain = self.args.chain

        if not Path("crystal_structure.pdb").exists():
            openawsem.helperFunctions.cleanPdb(
                [self.name],
                chain=chain,
                toFolder="cleaned_pdbs",
                verbose=self.args.verbose,
                keepIds=True,
                removeHeterogens=removeHeterogens
            )
            cleaned_pdb_path = Path(f"cleaned_pdbs/{self.pdb}")
            shutil.copy(cleaned_pdb_path,"crystal_structure.pdb")

        if chain == "-1":
            chain = openawsem.helperFunctions.getAllChains(
                "crystal_structure.pdb",
                removeDNAchains=True
            )
            print("Chains info read from crystal_structure.pdb, chains to simulate: ", chain)

        # for compute Q
        input_pdb_filename, cleaned_pdb_filename = openawsem.prepare_pdb(
            "crystal_structure.pdb",
            chain,
            use_cis_proline=False,
            keepIds=self.args.keepIds,
            removeHeterogens=removeHeterogens
        )
        openawsem.ensure_atom_order(input_pdb_filename)

        self.input_pdb_filename = input_pdb_filename
        self.cleaned_pdb_filename = cleaned_pdb_filename
        
        self.chain = openawsem.helperFunctions.getAllChains("crystal_structure-cleaned.pdb")
        openawsem.getSeqFromCleanPdb(input_pdb_filename, chains=chain, writeFastaFile=True)
        shutil.copy('crystal_structure.fasta',f'{self.name}.fasta')
        
        if self.args.extended:
            # print("Trying to create the extended structure extended.pdb using pymol, please ensure that pymol is installed and callable using 'pymol' in terminal.")
            # self.run_command(["python", f"{__location__}/helperFunctions/fasta2pdb.py", "extended", "-f", f"{self.name}.fasta"])
            # # print("If you has multiple chains, please use other methods to generate the extended structures.")
            # openawsem.helperFunctions.myFunctions.add_chain_to_pymol_pdb("extended.pdb")  # only work for one chain only now
            openawsem.helperFunctions.create_extended_pdb_from_fasta(f"{self.name}.fasta", output_file_name="extended.pdb")
            input_pdb_filename, cleaned_pdb_filename = openawsem.prepare_pdb("extended.pdb", "A", use_cis_proline=False, keepIds=args.keepIds, removeHeterogens=removeHeterogens)
            openawsem.ensure_atom_order(input_pdb_filename)
        
        shutil.copy('crystal_structure-cleaned.pdb',f'{self.pdb}')
        
        if self.args.keepLigands:
            # cleaned_pdb_filename = f"{name}-cleaned.pdb"
            # input_pdb_filename = f"{name}-openmmawsem.pdb"
            # do(f"grep 'ATOM' {input_pdb_filename} > tmp.pdb")
            # do(f"grep 'HETATM' {cleaned_pdb_filename} >> tmp.pdb")
            # do(f"mv tmp.pdb {input_pdb_filename}")

            shutil.copy("crystal_structure-cleaned.pdb", f"{self.name}-cleaned.pdb")
            with open("tmp.pdb", "w") as output_file:
                with open("crystal_structure-openmmawsem.pdb", "r") as input_file:
                    for line in input_file:
                        if "ATOM" in line:
                            output_file.write(line)
                with open("crystal_structure-cleaned.pdb", "r") as input_file:
                    for line in input_file:
                        if "HETATM" in line:
                            output_file.write(line)
            os.rename("tmp.pdb", f"{self.name}-openmmawsem.pdb")
        else:
            input_pdb_filename, cleaned_pdb_filename = openawsem.prepare_pdb(self.pdb, self.chain, keepIds=self.args.keepIds, removeHeterogens=removeHeterogens)
            openawsem.ensure_atom_order(input_pdb_filename)

    def generate_ssweight_from_stride(self):
        """
        Generate the secondary structure weight file (ssweight) using stride or Predict_Property.
        """
        self.run_command(["stride", "crystal_structure.pdb"], stdout="ssweight.stride")
        self.run_command(["python", __location__/"helperFunctions"/"stride2ssweight.py"], stdout="ssweight")
        protein_length = openawsem.helperFunctions.getFromTerminal("wc ssweight").split()[0]
        if int(protein_length) == 0:
            seq = openawsem.helperFunctions.read_fasta(f"{self.name}.fasta")
            protein_length = len(seq)
            print("impose no secondary bias.")
            print("you might want to install Predict_Property and use the predict_ssweight_from_fasta option.")
            with open("ssweight", "w") as out:
                for i in range(protein_length):
                    out.write("0.0 0.0\n")
        print(f"protein: {self.name}, length: {protein_length}")

    def generate_ssweight_from_fasta(self):

        # another option for secondary prediction bias generation is using "Predict_Property.sh -i {name}.fasta" to predict from fasta file.
        # but you need install it from https://github.com/realbigws/Predict_Property.
        # after installation, you can do the following to generate ssweight.
        # for me I put 'export Predict_Property="/Users/weilu/Research/Build/Predict_Property"' inside ~/.bash_profile file.
        self.run_command(["$Predict_Property/Predict_Property.sh", "-i", f"{self.name}.fasta"])
        from_secondary = f"{self.name}_PROP/{self.name}.ss3"
        toPre = "."
        to_ssweight = f"{toPre}/ssweight"
        print("convert ssweight")
        data = pd.read_csv(from_secondary, comment="#", names=["i", "Res", "ss3", "Helix", "Sheet", "Coil"], sep="\s+")
        with open(to_ssweight, "w") as out:
            for i, line in data.iterrows():
                if line["ss3"] == "H":
                    out.write("1.0 0.0\n")
                if line["ss3"] == "E":
                    out.write("0.0 1.0\n")
                if line["ss3"] == "C":
                    out.write("0.0 0.0\n")
                        
    def prepare_membrane_files(self):
        """
        Prepare the required membrane-related files (zim and zimPosition) if the membrane or hybrid options are specified.
        """
        self.run_command(["grep", "-E", "CB|CA  GLY", "crystal_structure-cleaned.pdb"], stdout="cbs.data")
        self.run_command(["awk", "{if($9>15) print \"1\"; else if($9<-15) print \"3\"; else print \"2\"}", "cbs.data"], stdout="zimPosition")

        openawsem.helperFunctions.create_zim(f"crystal_structure.fasta", tableLocation=__location__/"helperFunctions")

    def generate_fragment_memory(self,database = "cullpdb_pc80_res3.0_R1.0_d160504_chains29712",fasta = None,N_mem = 20,brain_damage = 1.0,fragmentLength = 9):
        """
        Generate the fragment memory file if the frag option is specified.
        """
        if fasta is None:
            fasta = f"{self.name}.fasta"

        #self.run_command(["cp", f"{__location__}/database/cullpdb_pc80_*", "."])
        # files_to_copy = [
        #     f"{database}",
        #     f"{database}.fasta",
        #     f"{database}.phr",
        #     f"{database}.pin",
        #     f"{database}.psq",
        # ]

        # for file_name in files_to_copy:
        #    shutil.copy(database/file_name, Path('.')/file_name)

        openawsem.helperFunctions.create_fragment_memories(database=database, fasta_file=fasta, memories_per_position=N_mem, 
                                                           brain_damage=brain_damage, fragment_length=fragmentLength, pdb_dir=openawsem.data_path.pdb, 
                                                           index_dir=openawsem.data_path.index, frag_lib_dir=openawsem.data_path.gro,
                                                           failed_pdb_list_file=openawsem.data_path.pdbfail, pdb_seqres=openawsem.data_path.pdbseqres,
                                                           weight=1, evalue_threshold=10000, cutoff_identical=90)

        # self.run_command([
        #     "python", __location__/"helperFunctions"/"MultCha_prepFrags_index.py",
        #     database, fasta, str(N_mem), str(brain_damage), str(fragmentLength)
        # ], stdout="logfile")

        # Check and correct the fragment memory file
        openawsem.helperFunctions.check_and_correct_fragment_memory("frags.mem")

        # Relocate the file to the fraglib folder
        openawsem.helperFunctions.relocate(fileLocation="frags.mem", toLocation="fraglib")

        # Replace the file path in frags.mem
        openawsem.helperFunctions.replace(f"frags.mem", f"{__location__}//Gros/", "./fraglib/")
        self.run_command(["cp", "frags.mem", "frag_memory.mem"])

    def generate_single_memory(self):
        for c in self.chain:
            # print(f"convert chain {c} of crystal structure to Gro file")
            self.run_command(["python", f"{__location__}/helperFunctions/Pdb2Gro.py", "crystal_structure-cleaned.pdb", f"{self.name}_{c}.gro", f"{c}"])
        
        seq_data = openawsem.helperFunctions.seq_length_from_pdb("crystal_structure-cleaned.pdb", self.chain)
        with open("single_frags.mem", "w") as out:
            out.write("[Target]\nquery\n\n[Memories]\n")
            for (chain_name, chain_start_residue_index, seq_length) in seq_data:
                # print(f"write chain {chain_name}")
                out.write(f"{self.name}_{chain_name}.gro {chain_start_residue_index} 1 {seq_length} 20\n")   # residue index in Gro always start at 1.

    def generate_charges(self):
        openawsem.helperFunctions.generate_charge_array(Path(f"{self.name}.fasta"),Path('charge.txt'))

    def copy_parameters(self, destination_folder='.'):
        # Copy the files using shutil.copy()
        shutil.copy(__location__/"parameters"/"burial_gamma.dat", destination_folder)
        shutil.copy(__location__/"parameters"/"gamma.dat", destination_folder)
        shutil.copy(__location__/"parameters"/"membrane_gamma.dat", destination_folder)
        shutil.copy(__location__/"parameters"/"anti_HB", destination_folder)
        shutil.copy(__location__/"parameters"/"anti_NHB", destination_folder)
        shutil.copy(__location__/"parameters"/"anti_one", destination_folder)
        shutil.copy(__location__/"parameters"/"para_HB", destination_folder)
        shutil.copy(__location__/"parameters"/"para_one", destination_folder)
    
    def copy_scripts(self,destination_folder='.'):
        """
        Copy required scripts and files to the current working directory.
        """
        mm_run_path = __location__ /"scripts"/ "mm_run.py"
        mm_analysis_path = __location__ /"scripts"/ "mm_analyze.py"
        forces_setup_path = __location__ /"scripts"/ "forces_setup.py"

        shutil.copy(mm_run_path, destination_folder)
        shutil.copy(mm_analysis_path, destination_folder)
        shutil.copy(forces_setup_path, destination_folder)


    def run(self):

        """
        Execute the main workflow of the AWSEMSimulationProject class.
        """
        
        #Log the command used to run this module.
        self.log_commandline_args()

        for protein in self.args.proteins:
            self.change_directory(self.base_folder)
            self.args.protein = protein
            project_folder = Path(os.path.splitext(os.path.basename(protein))[0])
            project_folder.mkdir(parents=True, exist_ok=True)
            
            # Prepare the input files
            if self.args.protein[-4:] == '.pdb':
                self.name, self.pdb = self.prepare_input_files_from_pdb(project_folder)
            elif self.args.protein[-6:] == ".fasta":
                self.name, self.pdb = self.prepare_input_files_from_fasta(project_folder)
            else:
                self.name, self.pdb = self.prepare_input_files_from_name(project_folder)
                
            print(self.name, self.pdb)
            
            #Change the directory
            with self.change_directory(project_folder):
            
                # Process the PDB files (clean, extract chains, generate extended structure)
                self.process_pdb_files()
                
                # Generate the secondary structure weight file (ssweight)
                if self.args.predict_ssweight_from_fasta:
                    self.generate_ssweight_from_fasta()
                else:
                    self.generate_ssweight_from_stride()
                
                # Prepare the membrane-related files if membrane or hybrid option is enabled
                if self.args.membrane or self.args.hybrid:
                    self.prepare_membrane_files()
                
                # Generate single memory file
                self.generate_single_memory()

                # Generate fragment memory files if the frag option is enabled
                if self.args.frag:
                    self.generate_fragment_memory(database=self.args.frag_database, fasta=self.args.frag_fasta, N_mem=self.args.frag_N_mem, brain_damage=self.args.frag_brain_damage, fragmentLength=self.args.frag_fragmentLength)

                #Generate charges
                self.generate_charges()
                
                # Copy required scripts to the current working directory
                self.copy_scripts()

                # Copy required parameters to the current working directory
                self.copy_parameters()

                print(f"{project_folder} project folder created")
                print("please modify the forces_setup.py if we want to change what energy terms to be used.")

import unittest
import tempfile
import shutil
from unittest.mock import MagicMock

class TestAWSEMSimulationProject(unittest.TestCase):

    def setUp(self):
        self.args = MagicMock()
        self.args.proteins = ['1r69']
        self.args.protein = '1r69'
        self.args.chain = "-1"
        self.args.debug = False
        self.args.keepLigands = False
        self.args.keepIds = False
        self.args.membrane = False
        self.args.hybrid = False
        self.args.frag = False
        self.args.predict_ssweight_from_fasta = False
        self.args.extended = False
        self.args.verbose = False
        
        self.temp_dir = Path(tempfile.mkdtemp())
        self.data_folder = __location__/'tests'/'data'
        self.project = AWSEMSimulationProject(self.args)
        self.project.base_folder = self.temp_dir



    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_prepare_input_files_from_pdb(self):
        pdb_path = Path(self.data_folder) / "1r69.pdb"
        self.args.protein = str(pdb_path)
        name, pdb = self.project.prepare_input_files_from_pdb(self.temp_dir)
        self.assertEqual(name, "1r69")
        self.assertEqual(pdb, "1r69.pdb")

    def test_prepare_input_files_from_fasta(self):
        fasta_path = Path(self.data_folder) / "1r69.fasta"
        self.args.protein = str(fasta_path)
        name, pdb = self.project.prepare_input_files_from_fasta(self.temp_dir)
        self.assertEqual(name, "1r69")
        self.assertEqual(pdb, "1r69.pdb")

    def test_prepare_input_files_from_name(self):
        self.args.protein = "1r69"
        name, pdb = self.project.prepare_input_files_from_name(self.temp_dir)
        self.assertEqual(name, "1r69")
        self.assertEqual(pdb, "1r69.pdb")

    def test_copy_scripts(self):
        self.project.copy_scripts(destination_folder=self.temp_dir)
        copied_files = ["mm_run.py", "mm_analysis.py", "forces_setup.py"]
        self.project.run_command(["ls"])
        for file in copied_files:
            print(self.project.base_folder/file)
            self.assertTrue((self.project.base_folder/file).exists())

    def test_run_command(self):
        with self.assertRaises(Exception):
            self.project.run_command(["ls", "non_existent_file"])
    

def main():
    args = parse_arguments()
    if args.test:
        unittest.main(argv=['first-arg-is-program-name'], exit=False) 
        exit()
    else:
        project = AWSEMSimulationProject(args)
        project.run()

if __name__=="__main__":
    main()