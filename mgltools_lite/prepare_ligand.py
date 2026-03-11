import os
from openbabel import openbabel

def prepare_ligand(input_file, output_pdbqt):
    ext = os.path.splitext(input_file)[1][1:]
    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats(ext, "pdbqt")
    mol = openbabel.OBMol()
    print(f"Procesando ligando: {input_file}")
    obc.ReadFile(mol, input_file)
    mol.AddHydrogens()
    builder = openbabel.OBBuilder()
    builder.Build(mol)
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(mol)
    obc.WriteFile(mol, output_pdbqt)
    print("✔ Ligando listo:", output_pdbqt)
