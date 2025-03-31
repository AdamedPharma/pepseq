import importlib
import rdkit
import rdkit.Chem


def mols_are_identical(
    mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol
) -> bool:
    """
    Check if two RDKit molecule objects are identical.
    This function compares two RDKit molecule objects to determine if they are
    identical by checking if each molecule is a substructure of the other,
    considering chirality.
    Args:
        mol1 (rdkit.Chem.rdchem.Mol): The first molecule to compare.
        mol2 (rdkit.Chem.rdchem.Mol): The second molecule to compare.
    Returns:
        bool: True if the molecules are identical, False otherwise.
    """
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def smiles_are_identical(smi1: str, smi2: str) -> bool:
    """
    Check if two SMILES strings represent the same molecule.
    Args:
        smi1 (str): The first SMILES string.
        smi2 (str): The second SMILES string.
    Returns:
        bool: True if the molecules are identical, False otherwise.
    """
    mol1 = rdkit.Chem.MolFromSmiles(smi1)
    mol2 = rdkit.Chem.MolFromSmiles(smi2)
    return mols_are_identical(mol1, mol2)


def smarts_are_identical(smi1: str, smi2: str) -> bool:
    """
    Check if two SMILES strings represent the same molecule.
    Args:
        smi1 (str): The first SMILES string.
        smi2 (str): The second SMILES string.
    Returns:
        bool: True if the molecules are identical, False otherwise.
    """
    mol1 = rdkit.Chem.MolFromSmarts(smi1)
    mol2 = rdkit.Chem.MolFromSmarts(smi2)
    if smi1 == smi2:
        return True
    else:
        return mols_are_identical(mol1, mol2)


def A_is_identical(A1, A2):
    all_keys = set(A1.keys()) | set(A2.keys())
    smiles_keys = ['smiles', 'smiles_radical']
    smarts_keys = ['smarts']
    other_keys = all_keys - set(smiles_keys)
    other_keys = other_keys - set(smarts_keys)

    for key in smiles_keys:
        if not smiles_are_identical(A1.get(key),A2.get(key)):
            return False
        
    for key in smarts_keys:
        if not smarts_are_identical(A1.get(key),A2.get(key)):
            return False
        
    for key in other_keys:
        if A1.get(key) != A2.get(key):
            return False
    return True


def aas_are_identical(aas1, aas2):
    all_keys = sorted( list( set(aas1.keys()) | set(aas2.keys()) ) )
    for key in all_keys:
        A1 = aas1.get(key)
        A2 = aas2.get(key)
        if not A_is_identical(A1, A2):
            return False
    return True


def check_db_json_some_keys(db_json1, db_json2):
    keys = ['coding', 'physico_chemical_properties', 'l_proteogenic',
           'd_proteogenic','aa_order', 'n_terms_order', 'c_terms_order',
           'protein_letters']
    for key in keys:
        if not db_json1.get(key) == db_json2.get(key):
            return False
    return True


def check_smi_j_other_keys(smi_j1, smi_j2):
    keys = ['n_terms', 'c_terms']
    for key in keys:
        if smi_j1.get(key) != smi_j2.get(key):
            return False
    return True


def check_db_json_all_keys(db_json1, db_json2):
    if not check_db_json_some_keys(db_json1, db_json2):
        return False
    
    smi_j = db_json1.get('smiles')
    smi_j_fx = db_json2.get('smiles')
    if not check_smi_j_other_keys(smi_j, smi_j_fx):
        return False
    smi_aa = smi_j.get('aa')
    smi_aa_fx = smi_j_fx.get('aa')
    
    if not aas_are_identical(smi_aa, smi_aa_fx):
        return False
    return True


def load_tests(name):
    """
    Load and yield tests from a specified module.
    Args:
        name (str): The name of the module to load tests from.
    Yields:
        object: Each test object from the module's `tests` variable.
    """
    # Load module which contains test data
    tests_module = importlib.import_module(name)
    # Tests are to be found in the variable `tests` of the module
    for test in tests_module.tests:
        yield test