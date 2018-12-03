# openmmawsem
An implementation of the AWSEM coarse-grained protein folding forcefield in OpenMM

Example:
1r69.(rename pdb file to crystal_structure.pdb)


steup simulation folder:
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

Additional function used can be found in https://github.com/luwei0917/awsemmd_script
