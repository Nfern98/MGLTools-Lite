from rdkit import Chem

AUTODOCK4_TYPES = {
    "C": "C","N":"N","O":"O","S":"S",
    "F":"F","CL":"Cl","BR":"Br","I":"I",
    "MG":"Mg","ZN":"Zn","CA":"Ca","FE":"Fe",
}

def autodock_atom_type(atom):
    sym=atom.GetSymbol().upper()
    if atom.GetIsAromatic():
        if sym=='C': return 'A'
        if sym=='N': return 'NA'
        if sym=='O': return 'OA'
        if sym=='S': return 'SA'
    return AUTODOCK4_TYPES.get(sym, sym)

def identify_hbond_roles(mol):
    donors=set(); acceptors=set()
    patt_d=Chem.MolFromSmarts('[N;H1,H2,+0,+1]')
    patt_a=Chem.MolFromSmarts('[$([O-]),$([O]=C),$([N]=C)]')
    for idx in range(mol.GetNumAtoms()):
        if mol.GetSubstructMatches(patt_d):
            for m in mol.GetSubstructMatches(patt_d):
                donors.add(m[0])
        if mol.GetSubstructMatches(patt_a):
            for m in mol.GetSubstructMatches(patt_a):
                acceptors.add(m[0])
    return donors,acceptors
