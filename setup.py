from setuptools import setup, find_packages

setup(
    name='mgltools-lite',
    version='1.3',
    packages=find_packages(),
    install_requires=['rdkit','openbabel-wheel'],
    description='Modern AutoDock4 ligand/receptor preparation tools',
)
