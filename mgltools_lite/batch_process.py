import os
from .prepare_ligand import prepare_ligand
from .prepare_receptor import prepare_receptor

def process_all_ligands(in_folder,out_folder):
    os.makedirs(out_folder, exist_ok=True)
    for f in os.listdir(in_folder):
        if f.lower().endswith(('pdb','mol','mol2')):
            prepare_ligand(os.path.join(in_folder,f),
                           os.path.join(out_folder, f.split('.')[0]+'.pdbqt'))

def process_all_receptors(in_folder,out_folder):
    os.makedirs(out_folder, exist_ok=True)
    for f in os.listdir(in_folder):
        if f.lower().endswith('pdb'):
            prepare_receptor(os.path.join(in_folder,f),
                             os.path.join(out_folder, f.split('.')[0]+'.pdbqt'))
