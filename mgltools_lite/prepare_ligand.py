from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from .autodock_typing import autodock_atom_type, identify_hbond_roles

def prepare_ligand(input_path, output_pdbqt):
    print(f"[Ligando] Procesando: {input_path}")

    # ------------------------------
    # Cargar con RDKit
    # ------------------------------
    mol = None
    ext = input_path.lower()

    if ext.endswith(".mol") or ext.endswith(".mol2"):
        mol = Chem.MolFromMolFile(input_path, removeHs=False)
    if mol is None and ext.endswith(".pdb"):
        mol = Chem.MolFromPDBFile(input_path, removeHs=False)

    if mol is None:
        raise ValueError(f"No se pudo cargar el ligando: {input_path}")

    # Añadir hidrógenos
    mol = Chem.AddHs(mol)

    # Generar conformación 3D
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)

    donors, acceptors = identify_hbond_roles(mol)

    # Guardar MOL temporal (necesario para OBMol)
    temp_mol = "temp_ligand.mol"
    Chem.MolToMolFile(mol, temp_mol)

    # ------------------------------
    # Convertir OpenBabel → PDBQT
    # ------------------------------
    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats("mol", "pdbqt")

    obmol = openbabel.OBMol()
    obc.ReadFile(obmol, temp_mol)

    # Cargas Gasteiger
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol)

    # ------------------------------
    # Aplicar TIPADO AD4 CORRECTO
    # ------------------------------
    for atom in openbabel.OBMolAtomIter(obmol):
        rd_atom = mol.GetAtomWithIdx(atom.GetIdx())
        atom.SetType(autodock_atom_type(rd_atom))

    # Torsiones
    obmol.SetTorsionEnergy(obmol.NumRotors())

    # Exportar PDBQT final
    obc.WriteFile(obmol, output_pdbqt)

    print(f"✔ Ligando preparado → {output_pdbqt}")
