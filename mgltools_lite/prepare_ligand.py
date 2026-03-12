from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from .autodock_typing import autodock_atom_type, identify_hbond_roles

def prepare_ligand(input_path, output_pdbqt):
    print(f"[Ligando] Procesando: {input_path}")

    # ----------- CARGAR RDKit -----------
    mol = None
    ext = input_path.lower()

    if ext.endswith(".mol") or ext.endswith(".mol2"):
        mol = Chem.MolFromMolFile(input_path, removeHs=False)
    elif ext.endswith(".pdb"):
        mol = Chem.MolFromPDBFile(input_path, removeHs=False)

    if mol is None:
        raise ValueError(f"No se pudo cargar el ligando: {input_path}")

    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)

    donors, acceptors = identify_hbond_roles(mol)

    # ----------- EXPORTAR RDKit → MOL TEMPORAL ------
    temp_mol = "temp_ligand.mol"
    Chem.MolToMolFile(mol, temp_mol)

    # ----------- CARGAR en OpenBabel ---------------
    obmol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("mol", "pdbqt")
    conv.ReadFile(obmol, temp_mol)

    # Cargas Gasteiger
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol)

    # ----------- TIPADO AD4 POR COORDENADAS ----------
    # RDKit atom mapping by position
    rd_coords = [(atom.GetIdx(), mol.GetConformer().GetAtomPosition(atom.GetIdx())) for atom in mol.GetAtoms()]

    for ob_atom in openbabel.OBMolAtomIter(obmol):

        # Buscar átomo RDKit más cercano por coordenadas
        ox, oy, oz = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()

        min_dist = 9999
        best_idx = None

        for idx, pos in rd_coords:
            d = (pos.x - ox)**2 + (pos.y - oy)**2 + (pos.z - oz)**2
            if d < min_dist:
                min_dist = d
                best_idx = idx

        rd_atom = mol.GetAtomWithIdx(best_idx)
        ad4_type = autodock_atom_type(rd_atom)

        ob_atom.SetType(ad4_type)

    # Torsiones AutoDock
    obmol.SetTorsionEnergy(obmol.NumRotors())

    conv.WriteFile(obmol, output_pdbqt)
    print(f"✔ Ligando PDBQT generado correctamente → {output_pdbqt}")
