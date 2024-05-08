# OpenAWSEM
## An implementation of the AWSEM coarse-grained protein folding forcefield in OpenMM



OpenAWSEM is an implementation of the AWSEM (Associative memory, Water-mediated Structure, and Energy Model) coarse-grained protein forcefield designed for use with the OpenMM simulation toolkit.

## Installation

### Conda

To install OpenAWSEM using Conda, execute the following command:

```bash
conda install -c conda-forge openawsem
```

### Git

This installation mode is recommended for users that want to contribute to the code and Wolynes lab members.

```bash
#Clone the awsem repository
git clone https://github.com/cabb99/openawsem.git
cd openawsem

# Create a new conda environment
conda create -n openawsem -c conda-forge --file requirements.txt
conda activate openawsem

# Install the package in editable mode
pip install -e .
```

## Requirements

### STRIDE
STRIDE is used for secondary structure prediction. 
Download and install STRIDE and add it to your PATH:
https://webclu.bio.wzw.tum.de/stride/
```bash
wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz
tar -xvzf stride.tar.gz
cd stride
make
echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc
```

### PSIBLAST    
Install psiblast using the distribution from bioconda:

```bash
conda install -c conda-forge -c bioconda blast
```

Alternatively Download and install psiblast and add it to your PATH: 
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

```bash
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/$(curl -s "https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/" | grep -o 'ncbi-blast-[0-9.]*+-x64-linux.tar.gz'| head -n 1)
tar -xvzf ncbi-*.tar.gz
cd ncbi*/bin
echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc
```

### PDB_SEQRES
* Download pdb_seqres.txt and put it in the cloned openawsem repository location

```bash
wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
OPENAWSEM_LOCATION=$(python -c "import openawsem; print(openawsem.__location__)")
cp pdb_seqres.txt $OPENAWSEM_LOCATION/data
```

### Predict_Property

For secondary structure prediction from the fasta file OpenAWSEM can use "Predict_Property.sh -i {name}.fasta".
Install it from https://github.com/realbigws/Predict_Property.
After installation, add Predict_property.sh to $PATH so it can be executed
For example add 'export PATH = $PATH:/Users/weilu/Research/Build/Predict_Property/' inside the ~/.bash_profile file.

## Configuration
OpenAWSEM allows users to configure data storage paths. To do this:

Create a .awsem directory in your home folder.
Inside .awsem, create a configuration file named config.ini to specify data paths. 
The default paths point to the local data directory inside the OpenAWSEM module.
Example config.ini:

```ini
[Data Paths]
blast = /home/USER/data/database/cullpdb_pc80_res3.0_R1.0_d160504_chains29712
gro = /home/USER/data/Gros
pdb = /home/USER/data/PDBs
index = /home/USER/data/Indices
pdbfail = /home/USER/data/notExistPDBsList
pdbseqres = /home/USER/data/pdb_seqres.txt
topology = /home/USER/topology
```

## Example
Simulation of the amino terminal domain of Phage 434 repressor (1r69)

1. **Activate the OpenMM Environment:**
   Activate the required environment for running simulations.
   ```bash
   source activate openmm
   ```

2. **Set Up the Simulation Folder:**
   Create a simulation folder using the `awsem_create` command. The awsem_create command will automatically download the corresponding pdb.
   ```bash
   awsem_create 1r69 --frag
   ```
   Alternatively, if you have the `1r69.pdb` file:
   ```bash
   awsem_create 1r69.pdb --frag
   ```

3. **Modify the forces_setup.py**

   The `forces_setup.py` script determines which force (energy) terms are included in the simulation. 
   To activate the fragment memory term uncomment the fragment memory term and comment the single memory term.
   ```python
      # templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),
   ```
   It should look like this:
   ```python
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=False),
      #  templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),
   ```
3. **Run the Simulation:**
   Execute the simulation using the `awsem_run` command, specifying the platform, number of steps, and start and end temperatures for the annealing simulation.
   As an example we are running 1e5 steps, but it is common to run from 5 to 30 million steps in a single run.
   
   ```bash
   awsem_run 1r69 --platform CPU --steps 1e5 --tempStart 800 --tempEnd 200 -f forces_setup.py
   ```

4. **Compute Energy and Q:**
   Analyze the simulation results and redirect the output to `info.dat`.
   ```bash
   awsem_analyze 1r69 > info.dat
   ```

5. **Run Local Scripts (Optional):**
   The scripts are copied to the project folder and can be modified as needed. To run the local scripts, use the following commands:
   ```bash
   ./mm_run.py 1r69 --platform CPU --steps 1e5 --tempStart 800 --tempEnd 200 -f forces_setup.py
   ./mm_analyze.py 1r69 > energy.dat
   ```

## Notes:
AWSEM is capable of modeling protein-DNA interactions when used together with open3SPN2, which can be found in a separate package at https://github.com/cabb99/open3spn2.

For small proteins, the LAMMPS version may be faster than OpenAWSEM, especially if a GPU is unavailable. Consider using http://awsem-md.org for such cases.

A quick check of the stability of a protein in AWSEM can be done using the frustratometer server http://frustratometer.qb.fcen.uba.ar/



## Data availability
Data related to the paper "OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations" is available at https://app.globus.org/file-manager?origin_id=b4cef8ce-7773-4016-8513-829f388f7986&origin_path=%2FopenAWSEM_data%2F

## Citation
Please cite the following paper when using OpenAWSEM:
Lu, W., Bueno, C., Schafer, N. P., Moller, J., Jin, S., Chen, X., ... & Wolynes, P. G. (2021). OpenAWSEM with Open3SPN2: A fast, flexible, and accessible framework for large-scale coarse-grained biomolecular simulations. PLoS computational biology, 17(2), e1008308. https://doi.org/10.1371/journal.pcbi.1008308
