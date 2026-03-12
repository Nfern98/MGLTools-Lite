import os
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import AllChem

def prepare_ligand(input_file, output_pdbqt):
    print(f"[Ligand] Procesando: {input_file}")

    # -------------------------------
    # 1) Leer con RDKit
    # -------------------------------
    ext = os.path.splitext(input_file)[1].lower()

    if ext == ".sdf":
        mol = Chem.SDMolSupplier(input_file, removeHs=False)[0]
    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(input_file, removeHs=False)
    elif ext == ".pdb":
        mol = Chem.MolFromPDBFile(input_file, removeHs=False)
    else:
        print("❌ Formato no reconocido, saltando:", input_file)
        return

    # -------------------------------
    # 2) Validación básica
    # -------------------------------
    if mol is None or mol.GetNumAtoms() < 2:
        print(f"❌ ERROR: El ligando '{input_file}' está vacío o roto (átomos={mol.GetNumAtoms() if mol else 0}). Saltando.")
        return

    # -------------------------------
    # 3) Generar conformación 3D si falta
    # -------------------------------
    if mol.GetNumConformers() == 0:
        print("   → Generando conformación 3D...")
        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol)
        except:
            print("❌ No se pudo generar conformación 3D, saltando.")
            return

    # -------------------------------
    # 4) Minimización segura
    # -------------------------------
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except:
        print("⚠ No se pudo minimizar, continuando sin minimización.")

    # -------------------------------
    # 5) Guardar como PDB temporal
    # -------------------------------
    tmp_pdb = output_pdbqt.replace(".pdbqt", ".tmp.pdb")
    Chem.MolToPDBFile(mol, tmp_pdb)

    # -------------------------------
    # 6) Convertir PDB → PDBQT con OpenBabel
    # -------------------------------
    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats("pdb", "pdbqt")

    obmol = openbabel.OBMol()
    obc.ReadFile(obmol, tmp_pdb)

    obmol.AddHydrogens()
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(obmol)

    obc.WriteFile(obmol, output_pdbqt)

    os.remove(tmp_pdb)

    print(f"✔ Ligando listo: {output_pdbqt}")
