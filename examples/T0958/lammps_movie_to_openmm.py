#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
from time import sleep
import fileinput

if(platform.system() == 'Darwin'):  # Mac system (local machine)
    OPENAWSEM_LOCATION = "/Users/mingchenchen/Documents/openmmawsem/openmmawsem/"
elif(platform.system() == 'Linux'):
    OPENAWSEM_LOCATION = '/projects/pw8/wl45/openmmawsem/'
else:
    print("system unknown")
sys.path.insert(0, OPENAWSEM_LOCATION)
from openmmawsem import *

a = open("end.pdb").read().split("END")
print(a[0])
import os
import fileinput
import platform
from openmmawsem import *
from myFunctions_helper import *

os.system("rm openmmMovie.pdb")
os.system(f"echo 'REMARK converted from awsem lammps output' >> openmmMovie.pdb")
for i in range(1):
    with open("tmp.pdb", "w") as out:
        out.write(a[i])
    input_pdb_filename, cleaned_pdb_filename =prepare_pdb("tmp.pdb", "T")
    os.system(f"echo 'MODEL  {i+1}' >> openmmMovie.pdb")
    os.system("cat tmp-openmmawsem.pdb >> openmmMovie.pdb")
