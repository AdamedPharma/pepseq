import json
import copy
import pkgutil
import importlib


from pepseq.tests_new.helpers import smiles_are_identical

from commands import (
    pepseq_to_smiles,
    calculate_json_from,
    read_smiles,
)

db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)

with open(db_path) as fp:
    db_json = json.load(fp)


def result_jsons_are_identical(j1: dict, j2: dict):
    """
    Compares two JSON-like dictionaries to determine if they are identical,
    excluding the "complete_smiles" key.
    Args:
        j1 (dict): The first JSON-like dictionary to compare.
        j2 (dict): The second JSON-like dictionary to compare.
    Returns:
        bool: True if the dictionaries are identical (excluding the "complete_smiles" key),
              and the values of the "complete_smiles" key are considered identical by
              the `smiles_are_identical` function. False otherwise.
    """
    j1_copy = copy.deepcopy(j1)
    j2_copy = copy.deepcopy(j2)

    smi1 = j1_copy.pop("complete_smiles")
    smi2 = j2_copy.pop("complete_smiles")

    if not smiles_are_identical(smi1, smi2):
        return False

    return j1_copy == j2_copy


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


# content of test_example.py
def pytest_generate_tests(metafunc):
    """
    A pytest hook to generate test cases dynamically based on fixture names.
    This function is called by pytest to generate test cases. It iterates over
    the fixture names provided in the test function and checks if any of them
    start with "data_". If a fixture name starts with "data_", it loads the
     associated test data using the `load_tests` function and then uses
     `metafunc.parametrize` to parameterize the test function with the loaded
     test data.
    Args:
        metafunc (Metafunc): The Metafunc object provided by pytest, which
                             contains information about the test function
                             being collected.
    """
    for fixture in metafunc.fixturenames:
        print(fixture)
        if fixture.startswith("data_"):
            # Load associated test data
            tests = load_tests(fixture)
            metafunc.parametrize(fixture, tests)


def test_peptide_from_pepseq_new(data_pepseq_smiles):
    """
    Test the conversion of a peptide sequence to SMILES format and validate the results.
    :param data_pepseq_smiles: A tuple containing the peptide sequence (pepseq)
         and the correct SMILES string (correct_smiles).
    :type data_pepseq_smiles: tuple

    :raises AssertionError: If the generated SMILES string does not match the
     correct SMILES string.

    The test performs the following steps:
    1. Converts the peptide sequence (pepseq) to a SMILES string.
    2. Reads the expected peptide sequences and modified SMILES strings from a file.
    3. Asserts that the read peptide sequences match the expected list.
    4. Splits the modified SMILES strings and asserts that each one matches the expected SMILES strings.
    5. Asserts that the generated SMILES string matches the correct SMILES string.
    """
    pepseq, correct_smiles = data_pepseq_smiles
    smiles = pepseq_to_smiles(pepseq)
    pepseq_list, mod_smiles_list = read_smiles("mypeptide.smi", "myppeptide_out")

    assert pepseq_list == ["H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"]
    mod_smiles_list = mod_smiles_list[0].split("\t")
    for i in range(len(mod_smiles_list)):
        assert smiles_are_identical(
            mod_smiles_list[i], ["[*:1]CNCC[*:2]", "[*:3]CNCCSP"][i]
        )

    assert smiles_are_identical(smiles, correct_smiles)


def test_calculate_json_from(data_calculate):
    """
    Test the `calculate_json_from` function to ensure it produces the expected JSON result.
    Args:
        data_calculate (tuple): A tuple containing the arguments and expected result.
            - args (list): The arguments to pass to the `calculate_json_from` function.
            - result (dict): The expected JSON result.
    Asserts:
        The JSON result from `calculate_json_from` matches the expected result.
    """
    args, result = data_calculate
    kwargs = {}
    assert result_jsons_are_identical(calculate_json_from(*args, **kwargs), result)
