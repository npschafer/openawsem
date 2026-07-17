"""An implementation of the AWSEM (Associative memory, Water-mediated Structure, and Energy Model) coarse-grained protein forcefield designed for use with the OpenMM simulation toolkit."""
from .openAWSEM import *
from pathlib import Path

__location__= Path(__file__).resolve().parent
from ._version import __version__

