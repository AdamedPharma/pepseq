import json
import os
from typing import Any, Dict

from pepseq.BuildingModifiedPeptideFromPeptideJSON import \
    get_smiles_from_peptide_json
from pepseq.BuildPeptideJSONFromSMILES import decompose_peptide_smiles
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.Peptide.models.Peptide import Peptide

absolute_path = os.path.dirname(__file__)
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def from_pepseq(pepseq: str, db_json: Dict = db_json) -> Peptide:
    """
    Convert a pepseq string to a Peptide object. Read peptide from PepSeq string.

    :param pepseq: The pepseq string to convert.
    :type pepseq: str
    :param db_json: The database JSON. Defaults to db_json.
    :type db_json: Dict

    :return: The converted Peptide object.
    :rtype: Peptide
    """

    mod_smiles = None
    peptide_json = get_pep_json(pepseq, db_json, mod_smiles)
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    peptide = Peptide(smiles, peptide_json)
    return peptide


def from_pepseq_and_mod_smiles(
    pepseq: str, mod_smiles: str, db_json: Dict = db_json
) -> Peptide:
    """
    Create a Peptide object from a peptide sequence and modified smiles.

    :param pepseq: The peptide sequence.
    :type pepseq: str
    :param mod_smiles: The modified smiles.
    :type mod_smiles: str
    :param db_json: The database JSON. Defaults to db_json.
    :type db_json: Dict

    :return: The created Peptide object.
    :rtype: Peptide
    """

    peptide_json = get_pep_json(pepseq, db_json, mod_smiles)
    peptide = from_json(peptide_json)
    return peptide


def from_smiles(smiles: str, db_json: Dict = db_json) -> Peptide:
    """
    Convert a SMILES string into a Peptide object.

    :param smiles: The SMILES string representing the peptide.
    :type smiles: str
    :param db_json: The database JSON containing information about the peptide. Defaults to db_json.
    :type db_json: Dict

    :return: The Peptide object representing the converted SMILES string.
    :rtype: Peptide
    """

    peptide_json = decompose_peptide_smiles(smiles, db_json)
    peptide = Peptide(smiles, peptide_json)
    return peptide


def from_json(peptide_json: Dict[str, Any], db: Dict = db_json) -> Peptide:
    """Read (modified) peptide from json

    pepseq json should looks as below:

    {"sequence":"H{Aib}EGTFTSDVSSYLEGQAAKEFIAWLVRGRG",
        "modifications":[
            {
            "modification_smiles": \
            "[*1]CC(=O)NCCCC(NC(C)=O)C(=O)Nc1ccc2oc(=O)cc(CC(=O)NCCOC" +
            "COCCC(=O)NCCCCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)c2c1[*2]",
            "connecting_residues": [17, 24]
            }
        ]
    }
    where [*1] and [*2] are attachement points

    :param peptide_json: json in format as above
    :type peptide_json: Dict[str, Any]
    :param residue_database: residue database. Defaults to db_json.
    :type residue_database: Dict

    :return: peptide object
    :rtype: Peptide

    Raises:
        Exception: _description_

    """

    smiles = get_smiles_from_peptide_json(peptide_json, db)
    peptide = Peptide(smiles, peptide_json)
    return peptide
