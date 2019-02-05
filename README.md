# openmmawsem
We are still in beta. all codes are subject to changes.

An implementation of the AWSEM coarse-grained protein folding forcefield in OpenMM

Installation
1. Download and install STRIDE and put it in your PATH: http://webclu.bio.wzw.tum.de/stride/
1. Download and install psiblast and put it in your PATH: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
1. Install conda or miniconda: https://conda.io/en/latest/miniconda.html
1. Create an openmm environment that includes the necessary Python packages.
```
conda create -n openmm python=3.6 biopython matplotlib pdbfixer mdtraj numpy pandas openmm
```
1. Set OPENAWSEM_LOCATION environment variable (the location where you cloned this Github repository).
```
export OPENAWSEM_LOCATION='/YOUR/OPENAWSEM/DIRECTORY/'
```
1. Activate the openmm environment.
```
source activate openmm
```

Example:
1r69.

setup simulation folder:
```
python3 YOUR_OPENAWSEM_LOCATION/mm_create_project.py 1r69 --frag
```

run the simulation:
```
python3 mm_run.py 1r69
```

compute energy and Q:
```
python3 mm_analysis.py 1r69 > energy.dat
```  
In comparison with the previous version, the q_value calculation inside the code has been modified in accordance with the default CACA Q as defined as Qw and Qo. There is a Qflag inside the newly added function (0 for Qw, 1 for Qo)
