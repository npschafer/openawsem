#!/usr/bin/env python2
# This script will take a STRIDE file and turn it into the standard format
# ssweight file for AWSEM-MD simulations.
# This is useful if you want to use the "correct" secondary structure
# assignments instead of a prediction.
# You can find stride files for all pdb files at:
# http://webclu.bio.wzw.tum.de/cgi-bin/stride/db.py
# The script assumes you have put the "plain" version of the stride
# database entry into a file called ssweight.stride.

# Nick Schafer 2/22/12

f = open('ssweight.stride', 'r')

empty_flag = True
for line in f:
    if not line.strip():
        continue
    else:
        empty_flag = False
        line=line.split()
        if line[0] == 'ASG':
            if line[6] == 'Strand':
                print('0.0 1.0')

            elif line[6] == 'AlphaHelix':
                print('1.0 0.0')

            else:
                print('0.0 0.0')

if empty_flag:
    raise IOError('ssweight.stride file was empty due to no hydrogen bonds in starting structure. Are you starting from extended structure? If so, use --predict_ssweight_from_fasta')
