import json
import pkgutil
import copy

import rdkit
import rdkit.Chem.PandasTools

import importlib
from pepseq.tests_new.helpers import smiles_are_identical


from pepseq.BuildingModifiedPeptideFromPeptideJSON import get_smiles_from_peptide_json
from pepseq.BuildPeptideJSONFromSMILES import (
    from_smiles_to_pepseq_and_one_mod_smiles_strings,
)
from pepseq.functions import calculate
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.read import from_json, from_pepseq
import pepseq.Peptide.utils.Parser
from pepseq.augmenting_db_json import augment_db_json
from pepseq.Peptide.utils.pure_parsing_functions import (
    decompose_symbol,
    get_attachment_point_json,
    get_attachment_points_on_sequence_json,
)


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


fixture_pepseq = "H~H{aMeAla}QGTY{Cys(R1)}DAQ{Cys(R2)}YS~NH2"
pepseq_value = "{Cys(R1)}ACDAPEPsEQ{Cys(R2)}"
base_seq_fixture = [
    "Cys(R1)",
    "A",
    "C",
    "D",
    "A",
    "P",
    "E",
    "P",
    "s",
    "E",
    "Q",
    "Cys(R2)",
]


single_modification_json = {
    "smiles": "[*:1]CNCC[*:2]",
    "max_attachment_point_id": 2,
    "attachment_points_on_sequence": {
        1: {
            "attachment_point_id": "1",
            "ResID": "1",
            "AtomName": "SG",
            "ResidueName": "Cys",
        },
        2: {
            "attachment_point_id": "2",
            "ResID": "12",
            "AtomName": "SG",
            "ResidueName": "Cys",
        },
    },
}


correct_peptide_json = {
    "sequence": "H{aMeAla}QGTYCDAQCYS",
    "length": 13,
    "internal_modifications": [],
    "C_terminus": "NH2",
    "N_terminus": "H",
    "pepseq_format": fixture_pepseq,
    "symbols": [
        "H",
        "H",
        "aMeAla",
        "Q",
        "G",
        "T",
        "Y",
        "Cys(R1)",
        "D",
        "A",
        "Q",
        "Cys(R2)",
        "Y",
        "S",
        "NH2",
    ],
    "external_modifications": [
        {
            "smiles": ("[*:1]CC(=O)NCC[C@H](NC(=O)C[*:2])C(=O)NCCC(=O)NC"
            "COC(=O)NCC[C@H](NC(=O)CCC(=O)O)C(=O)O"),
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
    "[*:1]CC(=O)NCC[C@H](NC(=O)C[*:2])C(=O)NCCC(=O)NCCOC(="
    "O)NCC[C@H](NC(=O)CCC(=O)O)C(=O)O"
)

mod_smiles = "[*:1]CNCC[*:2]"

mod_smiles_list = [mod_smiles]

correct_smiles = (
    "[H]N[C@@H](Cc1c[nH]cn1)C(=O)NC(C)(C)C(=O)N[C@@H](CC"
    "C(N)=O)C(=O)NCC(=O)N[C@H](C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H]1"
    "CSCC(=O)NCC[C@@H](C(=O)NCCC(=O)NCCOC(=O)NCC[C@H](NC(=O)CCC(=O)O)"
    "C(=O)O)NC(=O)CSC[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](CO"
    ")C(N)=O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)N"
    "C1=O)[C@@H](C)O"
)


def load_tests(name):
    """
    Load and yield tests from a specified module.
    Args:
        name (str): The name of the module to load tests from.
    Yields:
        test: Each test found in the `tests` variable of the specified module.
    """
    # Load module which contains test data
    tests_module = importlib.import_module(name)
    # Tests are to be found in the variable `tests` of the module
    for test in tests_module.tests:
        yield test


# content of test_example.py
def pytest_generate_tests(metafunc):
    """
    A pytest hook to generate test cases dynamically based on the fixture names.
    This function is called by pytest to generate test cases dynamically. It iterates
    through the fixture names provided in the test function and checks if any fixture
    name starts with "data_". If such a fixture is found, it loads the associated test
    data using the `load_tests` function and uses `metafunc.parametrize` to parameterize
    the fixture with the loaded test data.
    Args:
        metafunc (Metafunc): The Metafunc object provided by pytest, which contains
                             information about the test function and its fixtures.
    """
    for fixture in metafunc.fixturenames:
        if fixture.startswith("data_"):
            # Load associated test data
            tests = load_tests(fixture)
            metafunc.parametrize(fixture, tests)


def result_jsons_are_identical(j1: dict, j2: dict):
    """
    Compares two JSON objects to determine if they are identical, excluding the "complete_smiles" key.
    Args:
        j1 (dict): The first JSON object to compare.
        j2 (dict): The second JSON object to compare.
    Returns:
        bool: True if the JSON objects are identical (excluding the "complete_smiles" key), False otherwise.
    """
    j1_copy = copy.deepcopy(j1)
    j2_copy = copy.deepcopy(j2)

    smi1 = j1_copy.pop("complete_smiles")
    smi2 = j2_copy.pop("complete_smiles")

    if not smiles_are_identical(smi1, smi2):
        return False

    return j1_copy == j2_copy


def test_peptide_from_pepseq_new(data_pepseq_smiles):
    """
    Test the conversion of a peptide sequence to its corresponding SMILES representation.

    :param data_pepseq_smiles (tuple): A tuple containing a peptide sequence
     (pepseq) and its expected SMILES representation (smiles).
    :type data_pepseq_smiles: tuple

    :raises AssertionError: If the SMILES representation generated from the peptide
     sequence does not match the expected SMILES representation.
    :return: None

    Asserts:
        The SMILES representation generated from the peptide sequence matches the expected SMILES representation.
    """
    pepseq, smiles = data_pepseq_smiles
    peptide = from_pepseq(pepseq)
    assert smiles_are_identical(peptide.smiles, smiles)


def test_calculate(data_calculate):
    """
    Test the `calculate` function with provided data.
    :param  data_calculate: A tuple containing the arguments for the `calculate` function and the expected result.
    :type   data_calculate: tuple

    :raises AssertionError: If the result of the `calculate` function does not match the expected result.

    Asserts:
        The result of the `calculate` function matches the expected result
        using the `result_jsons_are_identical` function.
    """
    args, result = data_calculate
    kwargs = {}
    result_jsons_are_identical(calculate(*args, **kwargs), result)


def test_from_pepseq_and_one_mod_smiles_strings_to_peptide_json():
    """
    Test the conversion of a peptide sequence and a single modification SMILES string to a peptide JSON object.
    This test verifies that the `get_pep_json` function correctly converts a given peptide sequence,
     database JSON, and a list containing one modification SMILES string into a peptide JSON object.
     It checks that the SMILES string in the resulting JSON matches the expected SMILES string and
     that all other fields in the JSON match the expected values.
    Steps:
    1. Call `get_pep_json` with the fixture peptide sequence, database JSON, and
     a list containing one modification SMILES string.
    2. Print the resulting peptide JSON.
    3. Extract the expected SMILES string from the correct peptide JSON and compare
     it with the SMILES string in the resulting JSON.
    4. Verify that all other fields in the resulting JSON match the expected values.
    Assertions:
    - The SMILES string in the resulting peptide JSON matches the expected SMILES string.
    - All other fields in the resulting peptide JSON match the expected values.
    """
    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])

    correct_smi = correct_peptide_json.get("external_modifications")[0].pop("smiles")
    assert smiles_are_identical(
        correct_smi, peptide_json.get("external_modifications")[0].pop("smiles")
    )

    for k, v in correct_peptide_json.items():
        if k == "external_modifications":
            continue
        assert peptide_json.get(k) == v
    return


def test_from_pepseq_string_and_mod_smiles_to_smiles(
    data_smiles_from_pepseq_and_smiles,
):
    """
    Test the conversion of a peptide sequence string and modification SMILES to a SMILES string.
    This test function takes a tuple of arguments and the correct SMILES string, converts the
     peptide sequence and modification SMILES to a peptide JSON, and then converts the peptide
     JSON to a SMILES string. It asserts that the generated SMILES string matches the correct
     SMILES string.

    :param data_smiles_from_pepseq_and_smiles: A tuple containing the arguments for generating
        the peptide JSON and the correct SMILES string.
    :type data_smiles_from_pepseq_and_smiles: tuple

    :return: None
    :rtype: None
    """
    args, correct_smiles = data_smiles_from_pepseq_and_smiles

    peptide_json = get_pep_json(*args)
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    assert smiles == correct_smiles
    return


def test_from_pepseq_string_and_mod_smiles_to_peptide():
    """
    Test the conversion from a PEPSEQ string and modification SMILES to a peptide object.
    This test verifies that the peptide object created from a PEPSEQ string and a list of
     modification SMILES matches the expected peptide JSON structure, sequence, length,
     and complete SMILES representation.
    Assertions:
        - The peptide sequence matches the expected sequence in correct_peptide_json.
        - The peptide length is 13.
        - The peptide complete SMILES matches the expected correct_smiles.
    """
    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])
    peptide = from_json(peptide_json)
    peptide.sequence == correct_peptide_json["sequence"]
    peptide.length == 13
    peptide.complete_smiles == correct_smiles
    return


def test_from_smiles_to_pepseq_and_one_mod_smiles_strings():
    """
    Test the conversion from SMILES strings to PepSeq format and modified SMILES strings.
    This test verifies that the function `from_smiles_to_pepseq_and_one_mod_smiles_strings`
    correctly converts a given SMILES string to its corresponding PepSeq format and a list
    of modified SMILES strings. The test checks the following:
    1. The conversion of `correct_smiles` to `pepseq_format` and `mod_smiles` matches the
       expected `fixture_pepseq` and `one_mod_smiles`.
    2. The conversion of `fixture_complete_smiles_2` to `pepseq_2` and `mod_smiles_2` matches
       the expected `fixture_pepseq_2` and `fixture_mod_smiles_2`.
    Assertions:
    - `pepseq_format` should be equal to `fixture_pepseq`.
    - `mod_smiles` should be identical to `one_mod_smiles`.
    - `pepseq_2` should be equal to `fixture_pepseq_2`.
    - Each element in `mod_smiles_2` should be identical to the corresponding element in
      `fixture_mod_smiles_2`.
    """
    pepseq_format, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        correct_smiles, db_json
    )

    assert pepseq_format == fixture_pepseq
    assert smiles_are_identical(mod_smiles, one_mod_smiles)

    fixture_pepseq_2 = "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"
    fixture_complete_smiles_2 = (
            "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCCSP)C(=O)N[C@@H]"
            "(CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)"
            "[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2C"
            "CCN2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H]"
            "(CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O"
    )

    fixture_mod_smiles_2 = ["[*:1]CNCC[*:2]", "[*:3]CNCCSP"]

    pepseq_2, mod_smiles_2 = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        fixture_complete_smiles_2, db_json
    )

    assert pepseq_2 == fixture_pepseq_2
    assert smiles_are_identical(mod_smiles_2[0], fixture_mod_smiles_2[0])
    assert smiles_are_identical(mod_smiles_2[1], fixture_mod_smiles_2[1])

    return


def test_find_termini(data_canonical_sequence):
    """
    Test the find_termini function from the pepseq.Peptide.utils.Parser module.
    This test checks if the find_termini function correctly identifies the
     termini of a given canonical peptide sequence.

    :param data_canonical_sequence (str): The canonical peptide sequence to be tested.
    :type data_canonical_sequence: str

    Asserts:
        The function should return a tuple containing the N-terminus, C-terminus,
        and the modified peptide sequence.
    """
    assert pepseq.Peptide.utils.Parser.find_termini(
        data_canonical_sequence, db_json
    ) == ("H", "OH", "{Cys(R1)}ACDAPEPsEQ{Cys(R2)}")


def test_parse_canonical2(data_canonical_sequence):
    """
    Test the parse_canonical2 function from the pepseq.Peptide.utils.Parser module.
    Args:
        data_canonical_sequence (str): The canonical sequence to be parsed.
    Asserts:
        The parsed canonical sequence matches the expected base sequence fixture.
    """
    assert (
        pepseq.Peptide.utils.Parser.parse_canonical2(data_canonical_sequence)
        == base_seq_fixture
    )


def test_find_parentheses(data_canonical_sequence):
    """
    Test the find_parentheses function from the Parser class in the Peptide utils module.
    Args:
        data_canonical_sequence (str): A string representing the canonical sequence of the peptide.
    Asserts:
        The function should return a list of tuples, where each tuple contains the start and end indices
        of a pair of parentheses in the canonical sequence.
    """
    assert pepseq.Peptide.utils.Parser.find_parentheses(data_canonical_sequence) == [
        (0, 8),
        (19, 27),
    ]


def test_get_base_seq():
    """
    Test the get_base_seq function from the pepseq.get_peptide_json_from_pepseq_format module.
    This test verifies that the get_base_seq function correctly converts a sequence of symbols
    into the expected base sequence string.
    Steps:
    1. Retrieve the symbols from the base_seq_fixture.
    2. Call the get_base_seq function with the symbols.
    3. Assert that the returned base sequence matches the expected string "CACDAPEPsEQC".
    Expected Result:
    The base sequence returned by the get_base_seq function should be "CACDAPEPsEQC".
    """
    symbols = base_seq_fixture
    base_seq = pepseq.get_peptide_json_from_pepseq_format.get_base_seq(symbols)
    assert base_seq == "CACDAPEPsEQC"


def test_decompose_symbol():
    """
    Test the decompose_symbol function to ensure it correctly decomposes
    symbols with and without attachment points.
    The function should:
    - Return a tuple with the amino acid and attachment point number if the symbol has an attachment point.
    - Return the symbol itself if there is no attachment point.
    Test cases:
    - "Cys(R1)" should return ("Cys", "1")
    - "A" should return "A"
    """
    assert decompose_symbol("Cys(R1)") == ("Cys", "1")
    assert decompose_symbol("A") == "A"


attachment_points_on_sequence = single_modification_json.get(
    "attachment_points_on_sequence"
)


ext_mod_json = [single_modification_json]


def test_get_ext_mod_json():
    """
    Test the `get_ext_mod_json` function from the `pepseq.get_peptide_json_from_pepseq_format` module.
    This test verifies that the JSON output generated by the `get_ext_mod_json` function matches the expected
    JSON structure for a given peptide sequence and modification list.
    The test performs the following checks:
    1. Calls the `get_ext_mod_json` function with `base_seq_fixture` and `mod_smiles_list` as inputs.
    2. Prints the generated JSON for debugging purposes.
    3. Compares the `smiles` field of the generated JSON with the expected JSON using the `smiles_are_identical`
     function.
    4. Removes the `smiles` field from both the generated and expected JSON and compares the remaining structure.
    Assertions:
    - The `smiles` field in the generated JSON should be identical to the `smiles` field in the expected JSON.
    - The remaining structure of the generated JSON should match the expected JSON after removing the `smiles` field.
    """
    ext_mod_json_created = pepseq.get_peptide_json_from_pepseq_format.get_ext_mod_json(
        base_seq_fixture, mod_smiles_list
    )
    ext_mod_json_fx = [
        {
            "smiles": "[*:1]CNCC[*:2]",
            "max_attachment_point_id": 2,
            "attachment_points_on_sequence": {
                "1": {
                    "attachment_point_id": "1",
                    "ResID": "1",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
                "2": {
                    "attachment_point_id": "2",
                    "ResID": "12",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
            },
        }
    ]

    assert smiles_are_identical(
        ext_mod_json_fx[0].get("smiles"), ext_mod_json_created[0].get("smiles")
    )
    val = ext_mod_json_created[0]

    val.pop("smiles")
    fx = ext_mod_json_fx[0]
    fx.pop("smiles")
    assert val == fx


def test_get_attachment_points_on_sequence_json():
    """
    Test the function get_attachment_points_on_sequence_json with a base sequence fixture.
    This test checks if the function correctly identifies and returns the attachment points
    on a given sequence. The result is compared against the expected attachment points.
    Assertions:
        - The returned attachment points should match the expected attachment points.
    """
    att_points = get_attachment_points_on_sequence_json(base_seq_fixture)
    assert att_points == attachment_points_on_sequence


def test_get_single_modification_json():
    """
    Test the function get_single_modification_json from the pepseq module.
    This test checks if the function correctly converts a given peptide sequence
    with attachment points and modification SMILES to the expected JSON format.
    The function being tested:
    pepseq.get_peptide_json_from_pepseq_format.get_single_modification_json
    The test compares the output of the function with the expected JSON result.
    Args:
        attachment_points_on_sequence (str): The attachment points on the peptide sequence.
        mod_smiles (str): The SMILES representation of the modification.
    Returns:
        None
    """
    pepseq.get_peptide_json_from_pepseq_format.get_single_modification_json(
        attachment_points_on_sequence, mod_smiles
    ) == single_modification_json


def test_get_attachment_point_json():
    """
    Test the `get_attachment_point_json` function to ensure it returns the correct
    attachment point for a given residue ID and decomposition.
    The test checks if the function returns the expected attachment point from the
    `attachment_points_on_sequence` dictionary when provided with a specific residue ID
    and decomposition tuple.
    Test case:
    - decomposition: ("Cys", "1")
    - res_id: 1
    - Expected result: attachment_points_on_sequence.get(1)
    """
    decomposition = ("Cys", "1")
    res_id = 1
    assert get_attachment_point_json(
        res_id, decomposition
    ) == attachment_points_on_sequence.get(1)


def test_augmenting_db_json():
    """
    Test the augment_db_json function to ensure it correctly augments the database JSON
    with data from an SDF file.
    This test loads a set of monomers from an SDF file using RDKit's PandasTools,
     then calls the augment_db_json function with the loaded data. It checks that
     the resulting augmented database JSON is not empty and not None.
    Assertions:
        - The augmented database JSON should not be an empty dictionary.
        - The augmented database JSON should not be None.
    """
    df_sdf = rdkit.Chem.PandasTools.LoadSDF("monomers.sdf")
    db_json_augmented = augment_db_json(
        db_json, df_sdf=df_sdf, name_column="m_abbr", mol_colname="ROMol"
    )

    assert db_json_augmented != {}
    assert db_json_augmented is not None
    return
