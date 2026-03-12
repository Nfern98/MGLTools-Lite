# prepare_receptor.py

import os
from openbabel import openbabel

def prepare_receptor(input_pdb, output_pdbqt):

    print(f"[Receptor] Procesando: {input_pdb}")

    obc = openbabel.OBConversion()
    obc.SetInAndOutFormats("pdb", "pdbqt")

    mol = openbabel.OBMol()
    obc.ReadFile(mol, input_pdb)

    # -----------------------------------------------------------
    # 1. ELIMINAR AGUA (HOH / WAT)
    # -----------------------------------------------------------
    to_delete = []

    for atom in openbabel.OBMolAtomIter(mol):
        res = atom.GetResidue()
        if res and res.GetName().strip() in ["HOH", "WAT"]:
            to_delete.append(atom.GetId())

    for atom_id in to_delete:
        mol.DeleteAtom(mol.GetAtom(atom_id))

    # -----------------------------------------------------------
    # 2. AGREGAR HIDRÓGENOS POLARES
    # -----------------------------------------------------------
    mol.AddPolarHydrogens()

    # -----------------------------------------------------------
    # 3. CARGAS GASTEIGER
    # -----------------------------------------------------------
    charge_model = openbabel.OBChargeModel.FindType("gasteiger")
    charge_model.ComputeCharges(mol)

    # -----------------------------------------------------------
    # 4. ESCRIBIR ARCHIVO PDBQT
    # -----------------------------------------------------------
    obc.WriteFile(mol, output_pdbqt)

    print(f"✔ Receptor listo: {output_pdbqt}")
