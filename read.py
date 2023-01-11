import pkgutil
from typing import Any, Dict

from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.models.AminoAcidInstance import AminoAcidInstance
from Peptide.models.Peptide import Peptide
from Peptide.utils.Parser import parse_seq

db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
db = FileSystemDbRepo.read_from_json(db_path)


def from_pepseq(pepseq: str, residue_database: FileSystemDbRepo = db) -> Peptide:
    """Read peptide from PepSeq string"""
    symbols = parse_seq(pepseq, residue_database)
    n_terminus_code = symbols[0]
    c_terminus_code = symbols[-1]

    n_term_smiles = residue_database.n_terms_smi_codes[n_terminus_code]
    c_term_smiles = residue_database.c_terms_smi_codes[c_terminus_code]

    n_term_tuple = (n_terminus_code, n_term_smiles)
    c_term_tuple = (c_terminus_code, c_term_smiles)

    aa_symbols_list = symbols[1:-1]
    amino_acids = []

    for aa_symbol in aa_symbols_list:
        aai = AminoAcidInstance.MolFromSymbol(aa_symbol, residue_database)
        amino_acids.append(aai)

    amino_acids_tuple = (n_term_tuple, amino_acids, c_term_tuple)

    raise Exception("TO BE IMPLEMENTED ")


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
