"""
OpenAWSEM
Implementation of the AWSEM forcefield in openMM
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = "Calculates single residue frustration, and mutational frustration of proteins.".split("\n")[0]

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = None

setup(
    # Self-descriptive entries which should always be present
    name='openawsem',
    author='Carlos Bueno',
    author_email='carlos.bueno@rice.edu',
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,

    url='https://openawsem.org/',  # Website
    python_requires='>=3.8',
    install_requires=['biopython',
                      'matplotlib',
                      'numpy',
                      'pandas',
                      'openmm',
                      'pdbfixer',
                      'mdtraj',],              
    entry_points={
    'console_scripts': [
        'awsem = openawsem.scripts.awsem:main',
        'awsem_create = openawsem.scripts.mm_create_project:main',
        'awsem_run = openawsem.scripts.mm_run:main',
        'awsem_analyze = openawsem.scripts.mm_analyze:main',
    ],
    }
)
