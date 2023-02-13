import json
import pkgutil
from typing import Any, Dict

from BuildingModifiedPeptideFromPeptideJSON import get_smiles_from_peptide_json
from BuildPeptideJSONFromSMILES import decompose_peptide_smiles
from get_peptide_json_from_pepseq_format import get_pep_json
from Peptide.models.Peptide import Peptide

db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def from_pepseq(pepseq: str, db_json: Dict = db_json) -> Peptide:
    """Read peptide from PepSeq string"""
    mod_smiles = None
    peptide_json = get_pep_json(pepseq, db_json, mod_smiles)
    # peptide_json = get_pep_json(pepseq, mod_smiles)
    # peptide_json = get_peptide_json_from_sequence(pepseq, db_json)
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    peptide = Peptide(smiles, peptide_json)
    return peptide


def from_pepseq_and_mod_smiles(
    pepseq: str, mod_smiles: str, db_json: Dict = db_json
) -> Peptide:
    peptide_json = get_pep_json(pepseq, db_json, mod_smiles)
    peptide = from_json(peptide_json)
    return peptide


def from_smiles(smiles: str, db_json: Dict = db_json) -> Peptide:
    peptide_json = decompose_peptide_smiles(smiles, db_json)
    peptide = Peptide(smiles, peptide_json)
    return peptide


def from_json(peptide_json: Dict[str, Any], db_json: Dict = db_json) -> Peptide:
    """Read (modified) peptide from json

    pepseq json should looks as below:

    {"sequence":"H{Aib}EGTFTSDVSSYLEGQAAKEFIAWLVRGRG",
        "modifications":[
            {
            "modification_smiles": "[*1]CC(=O)NCCCC(NC(C)=O)C(=O)Nc1ccc2oc(=O)cc(CC(=O)NCCOCCOCCC(=O)NCCCCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)c2c1[*2]",
            "connecting_residues": [17, 24]
            }
        ]
    }
    where [*1] and [*2] are attachement points

    Args:
        json (Dict[str, Any]): json in format as above
        residue_database (FileSystemDbRepo, optional): residue database. Defaults to db.

    Raises:
        Exception: _description_

    Returns:
        Peptide: peptide object
    """
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    peptide = Peptide(smiles, peptide_json)
    return peptide
