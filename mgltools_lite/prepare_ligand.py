import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom
from rdkit.Chem.rdchem import Mol
from openbabel import openbabel

def prepare_ligand(input_file, output_pdbqt):
    print(f"[Ligando] Procesando: {input_file}")

    ext = os.path.splitext(input_file)[1].lower()

    # ------------------------------
    # 1. Cargar con RDKit (lectura robusta)
    # ------------------------------
    if ext == ".sdf":
        suppl = Chem.SDMolSupplier(input_file, removeHs=False)
        mol = suppl[0] if suppl and suppl[0] else None
    elif ext in [".mol2", ".mol"]:
        mol = Chem.MolFromMol2File(input_file, removeHs=False)
    elif ext == ".pdb":
        mol = Chem.MolFromPDBFile(input_file, removeHs=False)
    else:
        raise ValueError("Formato no soportado")

    if mol is None:
        print("⚠ RDKit no pudo cargar. Probando OpenBabel.")
        mol = None

    # ------------------------------
    # 2. Si RDKit falla → intentar OpenBabel
    # ------------------------------
    if mol is None:
        obc = openbabel.OBConversion()
        obc.SetInFormat(ext.replace(".", ""))
        obmol = openbabel.OBMol()
        obc.ReadFile(obmol, input_file)

        if obmol.NumAtoms() == 0:
            raise ValueError(f"❌ ERROR: {input_file} está vacío o corrupto.")

        # Convertir OBMol → RDKit
        obc.SetOutFormat("mol")
        temp_mol_file = "temp_ligand.mol"
        obc.WriteFile(obmol, temp_mol_file)
        mol = Chem.MolFromMolFile(temp_mol_file, removeHs=False)

        if mol is None:
            raise ValueError(f"❌ Fallo total al procesar {input_file}")

    # ------------------------------
    # 3. Reconstrucción de conectividad si falta
    # ------------------------------
    mol = Chem.AddHs(mol)

    # 3D embedding
    params = rdDistGeom.ETKDGv3()
    try:
        AllChem.EmbedMolecule(mol, params)
    except:
        print("⚠ Error embedding ETKDG, usando random coords")
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)

    # Minimización
    AllChem.MMFFOptimizeMolecule(mol)

    # ------------------------------
    # 4. Escribir MOL temporal y convertir a PDBQT con OpenBabel
    # ------------------------------
    temp_mol = "lig_temp.mol"
    Chem.MolToMolFile(mol, temp_mol)

    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats("mol", "pdbqt")
    obmol = openbabel.OBMol()
    obc.ReadFile(obmol, temp_mol)

    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol)

    obc.WriteFile(obmol, output_pdbqt)
    print(f"✔ Ligando exportado: {output_pdbqt}")utput_pdbqt}")
