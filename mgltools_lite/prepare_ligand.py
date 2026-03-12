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
    rd_coords = []
    conf = mol.GetConformer()

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        rd_coords.append((atom.GetIdx(), pos.x, pos.y, pos.z))

    for ob_atom in openbabel.OBMolAtomIter(obmol):
        ox, oy, oz = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()

        # Encontrar átomo RDKit más cercano por coordenadas
        best_idx = None
        best_dist = 1e9

        for idx, x, y, z in rd_coords:
            d = (x - ox)**2 + (y - oy)**2 + (z - oz)**2
            if d < best_dist:
                best_dist = d
                best_idx = idx

        rd_atom = mol.GetAtomWithIdx(best_idx)
        ad_type = autodock_atom_type(rd_atom)

        ob_atom.SetType(ad_type)

    # ----------- EXPORTAR LIGANDO FINAL ------------
    conv.WriteFile(obmol, output_pdbqt)

    print(f"✔ Ligando PDBQT generado correctamente → {output_pdbqt}")
