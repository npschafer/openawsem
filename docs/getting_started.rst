Getting Started
===============

Installation
------------

Conda
~~~~~
To install OpenAWSEM using Conda, execute the following command:

.. code-block:: bash

   conda install -c conda-forge openawsem

Git
~~~
This installation mode is recommended for users that want to contribute to the code and Wolynes lab members.

.. code-block:: bash

   # Clone the awsem repository
   git clone https://github.com/npschafer/openawsem.git
   cd openawsem

   # Create a new conda environment
   conda create -n openawsem -c conda-forge --file requirements.txt
   conda activate openawsem

   # Install the package in editable mode
   pip install -e .

Requirements
------------

STRIDE
~~~~~~
STRIDE is used for secondary structure prediction. 
Download and install STRIDE and add it to your PATH:
https://webclu.bio.wzw.tum.de/stride/

.. code-block:: bash

   mkdir stride
   cd stride
   wget https://webclu.bio.wzw.tum.de/stride/stride.tar.gz
   tar -xvzf stride.tar.gz
   make
   echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc

Note: If the webpage above becomes unavailable, please use an alternative repository like https://github.com/MDAnalysis/stride/tree/master/src .

PSIBLAST
~~~~~~~~
Install psiblast using the distribution from bioconda:

.. code-block:: bash

   conda install -c conda-forge -c bioconda blast

Alternatively Download and install psiblast and add it to your PATH: 
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

.. code-block:: bash

   wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/$(curl -s "https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/" | grep -o 'ncbi-blast-[0-9.]*+-x64-linux.tar.gz'| head -n 1)
   tar -xvzf ncbi-*.tar.gz
   cd ncbi*/bin
   echo 'export PATH=$PATH:'`pwd` >> ~/.bashrc

PDB_SEQRES
~~~~~~~~~~
* Download pdb_seqres.txt and put it in the cloned openawsem repository location

.. code-block:: bash

   wget https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt
   OPENAWSEM_LOCATION=$(python -c "import openawsem; print(openawsem.__location__)")
   cp pdb_seqres.txt $OPENAWSEM_LOCATION/data

Predict_Property
~~~~~~~~~~~~~~~~
For secondary structure prediction from the fasta file OpenAWSEM can use "Predict_Property.sh -i {name}.fasta".
Install it from https://github.com/realbigws/Predict_Property.
After installation, add Predict_property.sh to $PATH so it can be executed.
For example add 'export PATH=$PATH:/Users/weilu/Research/Build/Predict_Property/' inside the ~/.bash_profile file.


.. _example:

Example
-------

Simulation of the amino terminal domain of Phage 434 repressor (1r69)

1. **Activate the OpenMM Environment:**

   Activate the required environment for running simulations::

      source activate openmm

2. **Set Up the Simulation Folder:**

   Create a simulation folder using the ``awsem_create`` command. The command will automatically download the corresponding pdb::

      awsem_create 1r69 --frag

   Alternatively, if you have the ``1r69.pdb`` file::

      awsem_create 1r69.pdb --frag

3. **Modify the forces_setup.py**

   The ``forces_setup.py`` script determines which force (energy) terms are included in the simulation.
   To activate the fragment memory term, uncomment the fragment memory term and comment the single memory term.

   Original::

      # templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=True),
        templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),

   It should look like this::

        templateTerms.fragment_memory_term(oa, frag_file_list_file="./frags.mem", npy_frag_table="./frags.npy", UseSavedFragTable=False),
      #  templateTerms.fragment_memory_term(oa, frag_file_list_file="./single_frags.mem", npy_frag_table="./single_frags.npy", UseSavedFragTable=False),

4. **Run the Simulation:**

   Execute the simulation using the ``awsem_run`` command, specifying the platform, number of steps, and start and end temperatures for the annealing simulation.
   As an example, we are running 1e5 steps, but it is common to run from 5 to 30 million steps in a single run::

      awsem_run 1r69 --platform CPU --steps 1e5 --tempStart 800 --tempEnd 200 -f forces_setup.py

5. **Compute Energy and Q:**

   Analyze the simulation results and redirect the output to ``info.dat``::

      awsem_analyze 1r69 > info.dat

6. **Run Local Scripts (Optional):**

   The scripts are copied to the project folder and can be modified as needed. To run the local scripts, use the following commands::

      ./mm_run.py 1r69 --platform CPU --steps 1e5 --tempStart 800 --tempEnd 200 -f forces_setup.py
      ./mm_analyze.py 1r69 > energy.dat