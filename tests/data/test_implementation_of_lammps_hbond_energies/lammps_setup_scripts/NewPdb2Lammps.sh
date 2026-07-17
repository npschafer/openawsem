#!/bin/bash

# nifty one-liner from https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [[ ! $# -eq 2 ]]
then
	echo
	echo "> $0 pdb_id project_name"
	echo
	exit
fi

pdb_file=$1
output_file=$2

echo $pdb_file
echo $output_file

python2 ${SCRIPT_DIR}/PDBToCoordinates.py $pdb_file $pdb_file 
python2 ${SCRIPT_DIR}/CoordinatesToWorkLammpsDataFile.py $output_file".coord" "data."$output_file -b
