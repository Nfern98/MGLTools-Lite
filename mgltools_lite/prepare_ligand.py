import os
from openbabel import openbabel

def prepare_ligand(input_file, output_pdbqt):
    ext = os.path.splitext(input_file)[1][1:].lower()
    obc = openbabel.OBConversion()

    if ext in ["pdb", "mol2", "sdf"]:
        obc.SetInAndOutFormats(ext, "pdbqt")
    else:
        raise ValueError(f"Formato no soportado: {ext}")

    mol = openbabel.OBMol()
    print(f"[1/5] Cargando ligando: {input_file}")
    obc.ReadFile(mol, input_file)

    print(f"  → Átomos: {mol.NumAtoms()}")
    print(f"  → Conformers: {mol.NumConformers()}")

    # 2) Agregar hidrógenos
    print("[2/5] Agregando hidrógenos...")
    mol.AddHydrogens()

    # 3) Detectar si REALMENTE faltan coordenadas
    coords_missing = False
    for atom in openbabel.OBMolAtomIter(mol):
        x, y, z = atom.GetX(), atom.GetY(), atom.GetZ()
        if abs(x) < 1e-6 and abs(y) < 1e-6 and abs(z) < 1e-6:
            coords_missing = True
            break

    if coords_missing:
        print("[3/5] Faltan coordenadas → generando estructura 3D...")
        builder = openbabel.OBBuilder()
        builder.Build(mol)
    else:
        print("[3/5] Coordenadas correctas → NO se reconstruye el conformer")

    # 4) Cargas Gasteiger
    print("[4/5] Calculando cargas Gasteiger...")
    charge = openbabel.OBChargeModel.FindType("gasteiger")
    charge.ComputeCharges(mol)

    # 5) Exportar
    print("[5/5] Escribiendo PDBQT...")
    obc.WriteFile(mol, output_pdbqt)

    print(f"✔ Ligando preparado sin errores: {output_pdbqt}")
