from typing import TypeVar

from rdkit import Chem

from Peptide.models.Peptide import Peptide
from Peptide.utils.chemistry.MolExtended import MolObjExt
from Peptide.utils.Parser import read_sequence_txt, validate_sequence

Parser = TypeVar("Parser")
DataBase = TypeVar("DataBase")
SmilesParser = TypeVar("SmilesParser")
RepresentationFormat = TypeVar("RepresentationFormat")


class PeptideReader(object):
    def __init__(self):
        return


def read_sequence(
    sequence: str,
    db_api: DataBase,
) -> Peptide:
    valid = validate_sequence(sequence, db_api)

    try:
        n_term, amino_acids, c_term = read_sequence_txt(sequence, db_api)

    except:
        raise ()
    peptide = Peptide()
    peptide.amino_acids = (
        amino_acids  # setter can handle calculation of smiles and parameters
    )
    mol_obj_ext = MolObjExt(peptide.Mol)

    mol_obj_ext.n_term = n_term
    mol_obj_ext.c_term = c_term
    # consider  moving assigning n_term/c_term to Peptide.n_term.setter

    peptide.Mol = mol_obj_ext.mol

    peptide._n_term = n_term
    peptide._c_term = c_term

    return peptide


class SmilesReader(PeptideReader):
    def read(
        self,
        smiles: str,
        smiles_parser: SmilesParser,
        seq_parser: Parser,
        representation: RepresentationFormat = None,
        db_api=None,
    ) -> Peptide:
        mol = Chem.MolFromSmiles(smiles)
        seq = smiles_parser.read_smiles_txt(smiles, db_api)

        seq_reader = SeqReader()
        peptide = seq_reader.read(seq, seq_parser, db_api=db_api)

        c_mod_present = "OH"

        for c_mod in db_api.c_terms_order:
            seq_mod = "H~%s~%s" % (seq, c_mod)
            peptide_mod = seq_reader.read(seq_mod, seq_parser, db_api=db_api)
            is_mod = bool(mol.GetSubstructMatches(peptide_mod.Mol, useChirality=True))
            if is_mod:
                peptide = peptide_mod
                c_mod_present = c_mod

        for n_mod in db_api.n_terms_order:
            seq_mod = "%s~%s~%s" % (n_mod, seq, c_mod_present)
            peptide_mod = seq_reader.read(seq_mod, seq_parser, db_api=db_api)
            is_mod = bool(mol.GetSubstructMatches(peptide_mod.Mol, useChirality=True))
            if is_mod:
                peptide = peptide_mod
        return peptide
