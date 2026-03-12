import os
from .prepare_receptor import prepare_receptor
from .prepare_ligand import prepare_ligand

def process_all_ligands(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for file in os.listdir(input_folder):
        if file.lower().endswith(("pdb","mol2","sdf")):
            in_file = os.path.join(input_folder, file)
            out_file = os.path.join(output_folder, file.split('.')[0] + ".pdbqt")
            prepare_ligand(in_file, out_file)

def process_all_receptors(input_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    for file in os.listdir(input_folder):
        if file.lower().endswith("pdb"):
            in_file = os.path.join(input_folder, file)
            out_file = os.path.join(output_folder, file.split('.')[0] + ".pdbqt")
            prepare_receptor(in_file, out_file)
