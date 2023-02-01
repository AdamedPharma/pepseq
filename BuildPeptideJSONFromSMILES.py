import rdkit

from Backbone import BreakingIntoResidueCandidateSubgraphs, MarkingPeptideBackbone
from Peptide.utils.chemistry.ProcessResidueCandidateGraph import (
    decompose_residues_internal,
)


def decompose_peptide_smiles(smiles, cx_smarts_db):
    peptide_molecule = rdkit.Chem.MolFromSmiles(smiles)
    peptide_molecule = MarkingPeptideBackbone().execute(peptide_molecule)

    residues = BreakingIntoResidueCandidateSubgraphs().execute(peptide_molecule)

    seq, internal_modifications, external_modifications = decompose_residues_internal(
        residues, cx_smarts_db
    )
    return {
        "sequence": seq,
        "internal_modifications": internal_modifications,
        "external_modifications": external_modifications,
    }
