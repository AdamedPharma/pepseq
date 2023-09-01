import json
import pkgutil
import pandas as pd
import numpy as np
import networkx as nx

import rdkit
from pepseq.BuildingModifiedPeptideFromPeptideJSON import \
    BuildingModifiedPeptideFromPeptideJSON

from pepseq.Backbone import (BreakingIntoResidueCandidateSubgraphs,
                             MarkingPeptideBackbone)


from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


peptide_json = {
    "sequence": "CSCACGCK",
    "external_modifications": [
        {
            "smiles": "[1*]C([Br])CNP([Na])[2*]",
            "attachment_points_on_sequence": {
                1: {
                    "attachment_point_id": 1,
                    "ResID": "1",
                    "AtomName": "SG",
                    "ResidueName": "CYS",
                },
                2: {
                    "attachment_point_id": 2,
                    "ResID": "3",
                    "AtomName": "SG",
                    "ResidueName": "CYS",
                },
            },
        }
    ],
    "internal_modifications": {
        1: [
            {"ResID": 5, "AtomName": "SG", "ResidueName": "CYS"},
            {"ResID": 7, "AtomName": "SG", "ResidueName": "CYS"},
        ]
    },
}




def test_building():
    mol = BuildingModifiedPeptideFromPeptideJSON().execute(peptide_json, db_json)
    assert (
        rdkit.Chem.MolToSmiles(mol)
        == "[H]N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@@H](C)C(=O)N["
        + "C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(=O)O)NC(=O)CNC2=O)N"
        + "C(=O)[C@H](CO)NC1=O"
    )


def test_MarkingPeptideBackbone():
    mol_N_C_smiles_val = (
        "CC(=O)N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@"
        + "@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(N)=O)NC(=O)CNC"
        + "2=O)NC(=O)[C@H](CO)NC1=O"
    )

    peptide_molecule = rdkit.Chem.MolFromSmiles(mol_N_C_smiles_val)

    peptide_molecule2 = MarkingPeptideBackbone().execute(peptide_molecule)
    G = mol_to_nx(peptide_molecule2)
    selected_edges = [(u,v) for u,v,e in G.edges(data=True) if e.get('is_peptide_bond') == 'True']
    assert sorted(selected_edges) == [
        (16, 18), (21, 23), (30, 32), (42, 43), (46, 47), (49, 50), (55, 56)]


def test_BreakingIntoResidueCandidateSubgraphs():
    mol_N_C_smiles_val = (
            "CC(=O)N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@"
            + "@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(N)=O)NC(=O)CNC"
            + "2=O)NC(=O)[C@H](CO)NC1=O"
        )


    peptide_molecule = rdkit.Chem.MolFromSmiles(mol_N_C_smiles_val)

    subgraphs = BreakingIntoResidueCandidateSubgraphs().execute(peptide_molecule)

    ts =  tuple(
        [
            sum([ i[1]['atomic_num'] for i in list( subgraph.nodes(data=True) ) ]
               ) for subgraph in subgraphs])
    assert ts == (198, 33, 98, 65, 27, 41)