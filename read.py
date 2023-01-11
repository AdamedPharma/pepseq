import pkgutil

from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.models.Peptide import Peptide


db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
db = FileSystemDbRepo.read_from_json(db_path)


def from_pepseq(pepseq: str, residue_database: FileSystemDbRepo = db) -> Peptide:
    raise Exception("TO BE IMPLEMENTED ")


def from_smiles(smiles: str, residue_database: FileSystemDbRepo = db) -> Peptide:
    raise Exception("TO BE IMPLEMENTED ")
