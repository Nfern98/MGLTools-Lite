from setuptools import setup, find_packages

setup(
    name="mgltools-lite",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "openbabel-wheel",
        "rdkit-pypi"
    ],
    author="Nicolas Fernandez",
    description="Versión moderna de MGLTools Lite para preparar PDBQT.",
)
