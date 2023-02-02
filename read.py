import json
import pkgutil
from typing import Any, Dict

from BuildingModifiedPeptideFromPeptideJSON import (
    BuildingModifiedPeptideFromPeptideJSON,
    get_molecule_from_sequence,
)
from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.models.AminoAcidInstance import AminoAcidInstance
from Peptide.models.Peptide import Peptide
from Peptide.utils.Parser import parse_seq
from Peptide.utils.PeptideReader import read_sequence

db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
db = FileSystemDbRepo.read_from_json(db_path)
with open(db_path) as fp:
    db_json = json.load(db_path)


def from_pepseq(pepseq: str, db_json: Dict = db_json) -> Peptide:
    """Read peptide from PepSeq string"""
    mol = get_molecule_from_sequence(pepseq, db_json)
    return mol
    # return peptide


def from_smiles(smiles: str, residue_database: FileSystemDbRepo = db) -> Peptide:
    raise Exception("TO BE IMPLEMENTED ")


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
    return BuildingModifiedPeptideFromPeptideJSON.execute(peptide_json, db_json)
