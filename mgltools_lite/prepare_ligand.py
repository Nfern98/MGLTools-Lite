from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
from .autodock_typing import autodock_atom_type, identify_hbond_roles


def prepare_ligand(input_path, output_pdbqt):
    print(f"[Ligando] Procesando: {input_path}")

    # ============================================================
    # 1) Cargar el ligando SIEMPRE con OpenBabel primero
    #    → Esto repara ligandos rotos, sin enlaces, o con PDB malo
    # ============================================================
    obmol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInFormat("pdb")
    if not conv.ReadFile(obmol, input_path):
        raise ValueError(f"No se pudo leer el archivo PDB: {input_path}")

    # Reconstruir conectividad
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()

    # Guardar como MOL2 (incluye conectividad y tipos)
    conv.SetOutFormat("mol2")
    temp_mol2 = "temp_repaired.mol2"
    conv.WriteFile(obmol, temp_mol2)

    # ============================================================
    # 2) RDKit carga la molécula con conectividad correcta
    # ============================================================
    mol = Chem.MolFromMol2File(temp_mol2, removeHs=False)
    if mol is None:
        raise ValueError(f"No se pudo cargar en RDKit: {input_path}")

    # Añadir hidrógenos
    mol = Chem.AddHs(mol)

    # ============================================================
    # 3) EMBED MOLECULE + optimización segura
    # ============================================================
    res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if res != 0:
        print("⚠ Falló ETKDG, usando random embedding…")
        AllChem.EmbedMolecule(mol, randomSeed=0xBEEF)

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print("⚠ MMFF falló, optimizando con UFF…")
        AllChem.UFFOptimizeMolecule(mol)

    # ============================================================
    # 4) Guardar como MOL temporal para OpenBabel
    # ============================================================
    temp_mol = "temp_final.mol"
    Chem.MolToMolFile(mol, temp_mol)

    # ============================================================
    # 5) Abrir con OpenBabel y asignar cargas + tipado AD4
    # ============================================================
    conv.SetInAndOutFormats("mol", "pdbqt")
    obmol2 = openbabel.OBMol()
    conv.ReadFile(obmol2, temp_mol)

    # Cargas Gasteiger
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol2)

    # ============================================================
    # 6) Tipado AD4 por coordenadas (robusto)
    # ============================================================
    conf = mol.GetConformer()
    rd_coords = [
        (atom.GetIdx(),
         conf.GetAtomPosition(atom.GetIdx()).x,
         conf.GetAtomPosition(atom.GetIdx()).y,
         conf.GetAtomPosition(atom.GetIdx()).z)
        for atom in mol.GetAtoms()
    ]

    for ob_atom in openbabel.OBMolAtomIter(obmol2):
        ox, oy, oz = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()

        best_idx = None
        best_dist = 1e9

        # buscar átomo más cercano
        for idx, x, y, z in rd_coords:
            d = (x - ox)**2 + (y - oy)**2 + (z - oz)**2
            if d < best_dist:
                best_dist = d
                best_idx = idx

        rd_atom = mol.GetAtomWithIdx(best_idx)
        ad4 = autodock_atom_type(rd_atom)
        ob_atom.SetType(ad4)

    # ============================================================
    # 7) Guardar PDBQT final
    # ============================================================
    conv.WriteFile(obmol2, output_pdbqt)

    print(f"✔ Ligando preparado correctamente → {output_pdbqt}")
