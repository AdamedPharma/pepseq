import json
import pkgutil
import copy

import pytest
import rdkit
import rdkit.Chem.PandasTools


from pepseq.BuildingModifiedPeptideFromPeptideJSON import \
    get_smiles_from_peptide_json
from pepseq.BuildPeptideJSONFromSMILES import \
    from_smiles_to_pepseq_and_one_mod_smiles_strings
from pepseq.functions import calculate
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.read import from_json, from_pepseq
import pepseq.Peptide.utils.Parser
from pepseq.augmenting_db_json import augment_db_json
from pepseq.Peptide.utils.pure_parsing_functions import decompose_symbol, get_attachment_point_json, \
     get_attachment_points_on_sequence_json


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


fixture_pepseq = "H~H{aMeAla}QGTY{Cys(R1)}DAQ{Cys(R2)}YS~NH2"
pepseq_value = '{Cys(R1)}ACDAPEPsEQ{Cys(R2)}'
base_seq_fixture = ['Cys(R1)', 'A', 'C', 'D', 'A', 'P', 'E', 'P', 's', 'E', 'Q', 'Cys(R2)']


single_modification_json = {
        'smiles': '[1*]CNCC[2*]',
        'max_attachment_point_id': 2,
        'attachment_points_on_sequence': {
            1: {
                'attachment_point_id': '1',
                'ResID': '1',
                'AtomName': 'SG',
                'ResidueName': 'Cys'
            },
            2: {
                'attachment_point_id': '2',
                'ResID': '12',
                'AtomName': 'SG',
                'ResidueName': 'Cys'
            }
        }
    }

"""
single_modification_json = {
        'smiles': '[1*]CNCC[2*]',
        'max_attachment_point_id': 2,
        'attachment_points_on_sequence': {
            1: {
                'attachment_point_id': '1',
                'ResID': '1',
                'AtomName': 'SG',
                'ResidueName': 'Cys'
            },
            2: {
                'attachment_point_id': '2',
                'ResID': '12',
                'AtomName': 'SG',
                'ResidueName': 'Cys'
            }
        }
    }
"""

correct_peptide_json = {
    "sequence": "H{aMeAla}QGTYCDAQCYS",
    "length": 13,
    "internal_modifications": [],
    "C_terminus": "NH2",
    "N_terminus": "H",
    "pepseq_format": fixture_pepseq,
    "symbols": ["H", "H", "aMeAla", "Q", "G", "T", "Y", "Cys(R1)", "D", "A", "Q", "Cys(R2)", "Y", "S", "NH2"],
    "external_modifications": [
        {
            "smiles": "[1*]CC(=O)NCC[C@H](NC(=O)C[2*])C(=O)NCCC(=O)NC"
            + "COC(=O)NCC[C@H](NC(=O)CCC(=O)O)C(=O)O",
            "max_attachment_point_id": 2,
            "attachment_points_on_sequence": {
                1: {
                    "attachment_point_id": "1",
                    "ResID": "7",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
                2: {
                    "attachment_point_id": "2",
                    "ResID": "11",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
            },
        }
    ],
}


one_mod_smiles = (
    "[1*]CC(=O)NCC[C@H](NC(=O)C[2*])C(=O)NCCC(=O)NCCOC(="
    + "O)NCC[C@H](NC(=O)CCC(=O)O)C(=O)O"
)
mod_smiles = "[1*]CNCC[2*]"
mod_smiles_list = [mod_smiles]

correct_smiles = (
    "[H]N[C@@H](Cc1c[nH]cn1)C(=O)NC(C)(C)C(=O)N[C@@H](CC"
    + "C(N)=O)C(=O)NCC(=O)N[C@H](C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H]1"
    + "CSCC(=O)NCC[C@@H](C(=O)NCCC(=O)NCCOC(=O)NCC[C@H](NC(=O)CCC(=O)O)"
    + "C(=O)O)NC(=O)CSC[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](CO"
    + ")C(N)=O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)N"
    + "C1=O)[C@@H](C)O"
)






import importlib

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


def mols_are_identical(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol) -> bool:
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def smiles_are_identical(smi1: str, smi2: str) -> bool:
    mol1 = rdkit.Chem.MolFromSmiles(smi1)
    mol2 = rdkit.Chem.MolFromSmiles(smi2)
    return mols_are_identical(mol1, mol2)


def result_jsons_are_identical(j1: dict, j2: dict):
    j1_copy = copy.deepcopy(j1)
    j2_copy = copy.deepcopy(j2)

    smi1 = j1_copy.pop('complete_smiles')
    smi2 = j2_copy.pop('complete_smiles')

    if not smiles_are_identical(smi1, smi2):
        return False
    
    return j1_copy == j2_copy

def test_peptide_from_pepseq_new(data_pepseq_smiles):
    pepseq, smiles = data_pepseq_smiles
    peptide = from_pepseq(pepseq)
    assert smiles_are_identical(peptide.smiles, smiles)


def test_calculate(data_calculate):
    args, result = data_calculate
    result_jsons_are_identical(calculate(*args), result)


def test_from_pepseq_and_one_mod_smiles_strings_to_peptide_json():

    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])
    print(peptide_json)

    assert peptide_json == correct_peptide_json
    return


def test_from_pepseq_string_and_mod_smiles_to_smiles(data_smiles_from_pepseq_and_smiles):
    args, correct_smiles = data_smiles_from_pepseq_and_smiles

    peptide_json = get_pep_json( *args )
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    assert smiles == correct_smiles
    return


def test_from_pepseq_string_and_mod_smiles_to_peptide():
    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])
    peptide = from_json(peptide_json)
    peptide.sequence == correct_peptide_json["sequence"]
    peptide.length == 13
    peptide.complete_smiles == correct_smiles
    return


def test_from_smiles_to_pepseq_and_one_mod_smiles_strings():

    pepseq_format, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        correct_smiles, db_json
    )

    assert pepseq_format == fixture_pepseq
    assert mod_smiles == one_mod_smiles


    fixture_pepseq_2 = 'H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH'
    fixture_complete_smiles_2 = '[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCCSP)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O'
    fixture_mod_smiles_2 = ['[1*]CNCC[2*]', '[3*]CNCCSP']

    fixture_mod_smiles_3 = ['[*:1]CNCC[*:2]', '[*:3]CNCCSP']
    
    pepseq_2, mod_smiles_2 = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        fixture_complete_smiles_2, db_json
    )

    assert pepseq_2 == fixture_pepseq_2
    assert mod_smiles_2 == fixture_mod_smiles_2

    pepseq_3, mod_smiles_3 = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        fixture_complete_smiles_2, db_json, ketcher=True
    )

    assert pepseq_3 == fixture_pepseq_2
    assert mod_smiles_3 == fixture_mod_smiles_3

    return


def test_find_termini(data_canonical_sequence):
    assert pepseq.Peptide.utils.Parser.find_termini(data_canonical_sequence, db_json) == (
        'H', 'OH', '{Cys(R1)}ACDAPEPsEQ{Cys(R2)}')




def test_parse_canonical2(data_canonical_sequence):
    assert pepseq.Peptide.utils.Parser.parse_canonical2(data_canonical_sequence) == base_seq_fixture


def test_find_parentheses(data_canonical_sequence):
    assert pepseq.Peptide.utils.Parser.find_parentheses(data_canonical_sequence) == [
        (0, 8), (19, 27)]


def test_get_base_seq():
    symbols = base_seq_fixture
    base_seq = pepseq.get_peptide_json_from_pepseq_format.get_base_seq(symbols)
    assert base_seq == 'CACDAPEPsEQC'


def test_decompose_symbol():
    assert decompose_symbol('Cys(R1)') == ('Cys', '1')
    assert decompose_symbol('A') == 'A'




attachment_points_on_sequence = single_modification_json.get('attachment_points_on_sequence')


ext_mod_json = [ single_modification_json ]


def test_get_ext_mod_json():
    ext_mod_json_created = pepseq.get_peptide_json_from_pepseq_format.get_ext_mod_json(
    base_seq_fixture, mod_smiles_list)
    print(ext_mod_json_created)
    assert ext_mod_json_created == ext_mod_json



def test_get_attachment_points_on_sequence_json():
    att_points = get_attachment_points_on_sequence_json(base_seq_fixture)
    print(att_points)
    assert att_points == attachment_points_on_sequence


def test_get_single_modification_json():

    pepseq.get_peptide_json_from_pepseq_format.get_single_modification_json(
                    attachment_points_on_sequence, mod_smiles
                ) == single_modification_json


def test_get_attachment_point_json():
    decomposition = ('Cys', '1')
    res_id = 1
    assert get_attachment_point_json(res_id, decomposition
        ) == attachment_points_on_sequence.get(1)


def test_augmenting_db_json():
    df_sdf = rdkit.Chem.PandasTools.LoadSDF('monomers.sdf')
    db_json_augmented = augment_db_json(
        db_json, df_sdf=df_sdf, name_column = 'm_abbr', mol_colname='ROMol')

    assert bool(db_json_augmented) == True
    return

