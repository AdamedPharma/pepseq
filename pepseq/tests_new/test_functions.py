import json
import pkgutil
from pepseq.functions import calculate


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def test_calculate():
    pepseq_format = "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"
    smiles = ["[*:1]CNCC[*:2]", "[*:3]CNCC"]
    calculate(pepseq_format, smiles, db_json)
