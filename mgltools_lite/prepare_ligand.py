from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from .autodock_typing import autodock_atom_type


def prepare_ligand(input_path, output_pdbqt):
    print(f"[Ligando] Procesando: {input_path}")

    # =============================================================
    # 1) Cargar PDB con OpenBabel → eliminar enlaces duplicados
    # =============================================================
    obmol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInFormat("pdb")
    conv.ReadFile(obmol, input_path)

    # *** ESTA ES LA CLAVE ***
    # Reparar conectividad sin usar CONECT duplicados
    obmol.DeleteHydrogens()
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()

    # Guardar como MOL2
    conv.SetOutFormat("mol2")
    temp_mol2 = "temp_cleaned.mol2"
    conv.WriteFile(obmol, temp_mol2)

    # =============================================================
    # 2) RDKit carga el MOL2 limpio y válido
    # =============================================================
    mol = Chem.MolFromMol2File(temp_mol2, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit no pudo cargar el ligando: {input_path}")

    mol = Chem.AddHs(mol)

    # =============================================================
    # 3) Generar conformación 3D segura
    # =============================================================
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if result != 0:
        print("⚠ ETKDG falló, usando random embedding…")
        AllChem.EmbedMolecule(mol, randomSeed=123)

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        AllChem.UFFOptimizeMolecule(mol)

    # =============================================================
    # 4) Guardar como MOL para convertir a PDBQT
    # =============================================================
    temp_mol = "temp_ready.mol"
    Chem.MolToMolFile(mol, temp_mol)

    # Cargar en OpenBabel para tipado
    obmol2 = openbabel.OBMol()
    conv.SetInAndOutFormats("mol", "pdbqt")
    conv.ReadFile(obmol2, temp_mol)

    # Cargas Gasteiger
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol2)

    # =============================================================
    # 5) TIPADO AD4 por proximidad de coordenadas
    # =============================================================
    conf = mol.GetConformer()

    rd_coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        rd_coords.append((atom.GetIdx(), pos.x, pos.y, pos.z))

    for ob_atom in openbabel.OBMolAtomIter(obmol2):
        ox, oy, oz = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()

        best_idx = None
        best_dist = 1e9

        for idx, x, y, z in rd_coords:
            d = (x - ox)**2 + (y - oy)**2 + (z - oz)**2
            if d < best_dist:
                best_dist = d
                best_idx = idx

        rd_atom = mol.GetAtomWithIdx(best_idx)
        ob_atom.SetType(autodock_atom_type(rd_atom))

    # =============================================================
    # 6) Exportar como PDBQT
    # =============================================================
    conv.WriteFile(obmol2, output_pdbqt)

    print(f"✔ Ligando preparado correctamente → {output_pdbqt}")
