from openbabel import openbabel

def prepare_receptor(input_pdb, output_pdbqt):
    obc=openbabel.OBConversion()
    obc.SetInAndOutFormats("pdb","pdbqt")
    mol=openbabel.OBMol()
    obc.ReadFile(mol, input_pdb)

    for atom in list(mol.Atoms()):
        res=atom.GetResidue()
        if res and res.GetName().strip() in ['HOH','WAT']:
            res.DeleteAtom(atom)

    mol.AddPolarHydrogens()
    openbabel.OBChargeModel.FindType("gasteiger").ComputeCharges(mol)
    obc.WriteFile(mol, output_pdbqt)
