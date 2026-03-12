from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from .autodock_typing import autodock_atom_type, identify_hbond_roles

def prepare_ligand(input_path, output_pdbqt):
    print(f"[Ligando] Procesando: {input_path}")

    ext = input_path.lower()
    mol = None

    # -----------------------------
    # 1. Cargar MOL/MOL2 directamente
    # -----------------------------
    if ext.endswith(".mol") or ext.endswith(".mol2"):
        mol = Chem.MolFromMolFile(input_path, removeHs=False)

    # -----------------------------
    # 2. Si es PDB → RECONSTRUIR BONDS con OpenBabel
    # -----------------------------
    if mol is None and ext.endswith(".pdb"):
        obmol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat("pdb")
        conv.ReadFile(obmol, input_path)

        # Convertir a MOL2 (conectividad completa)
        conv.SetOutFormat("mol2")
        temp_mol2 = "temp_fixed.mol2"
        conv.WriteFile(obmol, temp_mol2)

        mol = Chem.MolFromMol2File(temp_mol2, removeHs=False)

    if mol is None:
        raise ValueError(f"No se pudo cargar el ligando {input_path}")

    # -----------------------------
    # 3. Añadir hidrógenos
    # -----------------------------
    mol = Chem.AddHs(mol)

    # -----------------------------
    # 4. Embedding seguro 3D
    # -----------------------------
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        print("⚠ Falló el embedding, intentando random...")
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)

    AllChem.MMFFOptimizeMolecule(mol)

    # -----------------------------
    # 5. Escribir MOL temporal
    # -----------------------------
    temp_mol = "temp_clean.mol"
    Chem.MolToMolFile(mol, temp_mol)

    # -----------------------------
    # 6. Convertir a OBMol para PDBQT
    # -----------------------------
    conv.SetInAndOutFormats("mol", "pdbqt")
    obmol = openbabel.OBMol()
    conv.ReadFile(obmol, temp_mol)

    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol)

    # -----------------------------
    # 7. Tipado AD4 por coordenadas
    # -----------------------------
    rd_coords = []
    conf = mol.GetConformer()

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        rd_coords.append((atom.GetIdx(), pos.x, pos.y, pos.z))

    for ob_atom in openbabel.OBMolAtomIter(obmol):
        ox, oy, oz = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()

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

    # -----------------------------
    # 8. Guardar PDBQT final
    # -----------------------------
    conv.WriteFile(obmol, output_pdbqt)

    print(f"✔ Ligando listo → {output_pdbqt}")
