from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from .autodock_typing import autodock_atom_type, identify_hbond_roles

def prepare_ligand(input_path, output_pdbqt):
    mol=None
    if input_path.lower().endswith('.mol') or input_path.lower().endswith('.mol2'):
        mol=Chem.MolFromMolFile(input_path, removeHs=False)
    if mol is None and input_path.lower().endswith('.pdb'):
        mol=Chem.MolFromPDBFile(input_path, removeHs=False)
    if mol is None:
        raise ValueError("Cannot read ligand file")

    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    donors,acceptors=identify_hbond_roles(mol)

    temp="temp_ligand.mol"
    Chem.MolToMolFile(mol, temp)

    obc=openbabel.OBConversion()
    obc.SetInAndOutFormats("mol","pdbqt")
    obmol=openbabel.OBMol()
    obc.ReadFile(obmol, temp)

    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol)

    for atom in obmol.Atoms():
        rd_atom = mol.GetAtomWithIdx(atom.GetIdx())
        atom.SetType(autodock_atom_type(rd_atom))

    obmol.SetTorsionEnergy(obmol.NumRotors())
    obc.WriteFile(obmol, output_pdbqt)
