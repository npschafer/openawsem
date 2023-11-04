#!/bin/bash

env_name="test_env"
setup_py_path="setup.py"
requirements_path="requirements.txt"
env_path="$HOME/miniconda3/envs/$env_name"

echo "Checking current environment..."
if [[ $CONDA_DEFAULT_ENV == $env_name ]]; then
    echo "Deactivating current environment..."
    conda deactivate
fi

echo "Removing environment '$env_name' if it exists..."
conda env remove --name $env_name -y

echo "Creating new environment '$env_name'..."
conda create --prefix $env_path --file $requirements_path -y

echo "Activating environment '$env_name'..."
source activate $env_name

echo "Checking Python path..."
which python

echo "Installing using setup.py..."
pip install -e .

echo "Testing completed."

