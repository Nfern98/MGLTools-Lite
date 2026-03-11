import os
from openbabel import openbabel

def prepare_receptor(input_pdb, output_pdbqt):
    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats("pdb", "pdbqt")
    mol = openbabel.OBMol()
    print(f"Procesando receptor: {input_pdb}")
    obc.ReadFile(mol, input_pdb)

    for atom in openbabel.OBMolAtomIter(mol):
        res = atom.GetResidue()
        if res and res.GetName().strip() in ["HOH","WAT"]:
            res.DeleteAtom(atom)

    mol.AddPolarHydrogens()
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(mol)
    obc.WriteFile(mol, output_pdbqt)
    print("✔ Receptor listo:", output_pdbqt)
