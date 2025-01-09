import json
import pkgutil
from pepseq.functions import calculate


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def test_calculate():
    """
    Test the `calculate` function with a given peptide sequence format and SMILES strings.
    This test checks if the `calculate` function correctly processes the provided
    peptide sequence format and SMILES strings, and interacts with the database
    JSON object.
    The peptide sequence format used in this test is:
    "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"
    The SMILES strings used in this test are:
    ["[*:1]CNCC[*:2]", "[*:3]CNCC"]
    Args:
        None
    Returns:
        None
    """
    pepseq_format = "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"
    smiles = ["[*:1]CNCC[*:2]", "[*:3]CNCC"]
    calculate(pepseq_format, smiles, db_json)
