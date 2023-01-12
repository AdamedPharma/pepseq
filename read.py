import pkgutil
from typing import Any, Dict

from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.models.AminoAcidInstance import AminoAcidInstance
from Peptide.models.Peptide import Peptide
from Peptide.utils.Parser import parse_seq
from Peptide.utils.PeptideReader import read_sequence

db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
db = FileSystemDbRepo.read_from_json(db_path)


def from_pepseq(pepseq: str, residue_database: FileSystemDbRepo = db) -> Peptide:
    """Read peptide from PepSeq string"""
    peptide = read_sequence(pepseq, residue_database)
    return peptide


def from_smiles(smiles: str, residue_database: FileSystemDbRepo = db) -> Peptide:
    raise Exception("TO BE IMPLEMENTED ")


def from_json(json: Dict[str, Any], residue_database: FileSystemDbRepo = db) -> Peptide:
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
    raise Exception("TO BE IMPLEMENTED")
