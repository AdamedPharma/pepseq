import json, copy
import pkgutil
import rdkit
import importlib


from pepseq.tests_new.helpers import mols_are_identical, smiles_are_identical

from commands import pepseq_to_smiles, calculate_json_from, read_smiles, \
    augment_db_json_command

db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)

with open(db_path) as fp:
    db_json = json.load(fp)

# from_smiles_to_pepseq_and_one_mod_smiles_strings

def result_jsons_are_identical(j1: dict, j2: dict):
    j1_copy = copy.deepcopy(j1)
    j2_copy = copy.deepcopy(j2)

    smi1 = j1_copy.pop('complete_smiles')
    smi2 = j2_copy.pop('complete_smiles')

    if not smiles_are_identical(smi1, smi2):
        return False
    
    return j1_copy == j2_copy


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
    #pepseq_list, # read_smiles('mypeptide.smi', 'myppeptide_out')
    # augment_db_json_command('my_monomers.sdf', 'augmented_db.json')    
    #peptide = from_pepseq(pepseq)
    pepseq_list, mod_smiles_list = read_smiles('mypeptide.smi', 'myppeptide_out')

    assert pepseq_list == ['H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH']
    assert mod_smiles_list == ['[1*]CNCC[2*]\t[3*]CNCCSP']

    assert smiles_are_identical(smiles, correct_smiles)


def test_calculate_json_from(data_calculate):
    args, result = data_calculate
    #result_jsons_are_identical(calculate(*args), result)
    # calculate_json_from('CH3~SC{R1}AFC~NH2', '[*1]CCC')
    kwargs = {}
    kwargs['ketcher'] = True
    assert result_jsons_are_identical(calculate_json_from(*args, **kwargs), result)

"""
def test_calculate_json_from_old_format(data_calculate_old_format):
    args, result = data_calculate
    #result_jsons_are_identical(calculate(*args), result)
    # calculate_json_from('CH3~SC{R1}AFC~NH2', '[*1]CCC')
    kwargs = {}
    kwargs['ketcher'] = True
    assert result_jsons_are_identical(calculate_json_from(*args, **kwargs), result)
"""


def test_ketchers():
    res_ketcher_true = calculate_json_from('CH3~SC{Cys(R1)}AFC~NH2', '[*:1]CCC', ketcher=True)
    res_ketcher_false = calculate_json_from('CH3~SC{Cys(R1)}AFC~NH2', '[1*]CCC', ketcher=False)

    complete_smiles_fx = 'CCCSC[C@H](NC(=O)[C@H](CS)NC(=O)[C@H](CO)NC)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CS)C(N)=O'

    assert smiles_are_identical(res_ketcher_true.get('complete_smiles'), complete_smiles_fx)
    assert smiles_are_identical(res_ketcher_false.get('complete_smiles'), complete_smiles_fx)

    assert res_ketcher_true.get('length') == 6
    assert res_ketcher_false.get('length') == 6

    assert res_ketcher_true.get('mw') == 687.91
    assert res_ketcher_false.get('mw') == 687.91

    assert res_ketcher_true.get('sequence') == 'SCCAFC'
    assert res_ketcher_false.get('sequence') == 'SCCAFC'
    
    
    res_ketcher_true_nd = calculate_json_from('CH3~SC{Cys(R1)}AFC~NH2', '[*:1]CCC')
    res_ketcher_false_nd = calculate_json_from('CH3~SC{Cys(R1)}AFC~NH2', '[1*]CCC')

    
    assert smiles_are_identical(res_ketcher_true_nd.get('complete_smiles'), complete_smiles_fx)
    assert smiles_are_identical(res_ketcher_false_nd.get('complete_smiles'), complete_smiles_fx)

    assert res_ketcher_true_nd.get('length') == 6
    assert res_ketcher_false_nd.get('length') == 6

    assert res_ketcher_true_nd.get('mw') == 687.91
    assert res_ketcher_false_nd.get('mw') == 687.91

    assert res_ketcher_true_nd.get('sequence') == 'SCCAFC'
    assert res_ketcher_false_nd.get('sequence') == 'SCCAFC'
    
    return


