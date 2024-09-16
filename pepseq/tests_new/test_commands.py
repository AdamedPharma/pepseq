import json
import pkgutil
import rdkit

from commands import pepseq_to_smiles, calculate_json_from, read_smiles, \
    augment_db_json_command

db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)

with open(db_path) as fp:
    db_json = json.load(fp)


import importlib

def mols_are_identical(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol) -> bool:
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def smiles_are_identical(smi1: str, smi2: str) -> bool:
    mol1 = rdkit.Chem.MolFromSmiles(smi1)
    mol2 = rdkit.Chem.MolFromSmiles(smi2)
    return mols_are_identical(mol1, mol2)


def load_tests(name):
    # Load module which contains test data
    tests_module = importlib.import_module(name)
    # Tests are to be found in the variable `tests` of the module
    for test in tests_module.tests:
        yield test


# content of test_example.py
def pytest_generate_tests(metafunc):

    for fixture in metafunc.fixturenames:
        print(fixture)
        if fixture.startswith('data_'):
            # Load associated test data
            tests = load_tests(fixture)
            metafunc.parametrize(fixture, tests)







def test_peptide_from_pepseq_new(data_pepseq_smiles):
    pepseq, correct_smiles = data_pepseq_smiles
    smiles = pepseq_to_smiles(pepseq)
    # calculate_json_from('CH3~SC{R1}AFC~NH2', '[*1]CCC')
    # read_smiles('mypeptide.smi', 'myppeptide_out')
    # augment_db_json_command('my_monomers.sdf', 'augmented_db.json')    
    #peptide = from_pepseq(pepseq)
    assert smiles_are_identical(smiles, correct_smiles)