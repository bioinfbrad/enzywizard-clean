#!/usr/bin/env python
from setuptools import setup, find_packages
import os

# Read the version from version.py without importing the package
version_file = os.path.join(os.path.dirname(__file__), 'src', 'enzywizard_clean', 'version.py')
with open(version_file) as f:
    exec(f.read())  # defines __version__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzywizard-clean",
    version=__version__,                     # dynamically read from version.py (1.0.1)
    author="bioinfbrad",
    description=(
        "Clean protein structures, generate multi-format protein files (CIF, PDB, and FASTA), "
        "and provide a detailed traceable cleaning report."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioinfbrad/enzywizard-clean",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.86",          # For structure I/O, residue handling
        "openmm>=8.0",              # Molecular mechanics engine (used for hydrogen addition)
        "pdbfixer>=1.12",           # PDBFixer APIs for cleaning
        "numpy>=1.23.5",            # Numerical operations
        "scipy>=1.15.2",            # Scientific calculations (used by other dependencies)
        "packaging",                # Version handling (used internally)
    ],
    entry_points={
        "console_scripts": [
            "enzywizard-clean = enzywizard_clean.cli:main",
        ],
    },
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
