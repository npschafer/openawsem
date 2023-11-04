from .openAWSEM import *
from . import _version
from pathlib import Path

__location__= Path(__file__).resolve().parent
__version__ = _version.get_versions()['version']

