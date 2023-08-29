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


from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol


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


def json_to_nx(mol_json: dict):
    """
    transforms custom JSON format
    into nx.classes.graph.Graph graph
    """

    G = nx.Graph()

    chirality_encoding = {
        0: rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
        1: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
        2: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    }

    hybridization_encoding = {
        1: rdkit.Chem.rdchem.HybridizationType.S,
        2: rdkit.Chem.rdchem.HybridizationType.SP,
        3: rdkit.Chem.rdchem.HybridizationType.SP2,
        4: rdkit.Chem.rdchem.HybridizationType.SP3,
    }

    bond_type_encoding = {
        1: rdkit.Chem.rdchem.BondType.SINGLE,
        2: rdkit.Chem.rdchem.BondType.DOUBLE
    }

    for (atomic_num, formal_charge, chiral_tag_encoded, hybridization_encoded,
         num_explicit_hs, is_aromatic, isotope, AtomName,
         ResID, node_id) in nodes_tuple:
        G.add_node(
            node_id,
            atomic_num=atomic_num,
            formal_charge=formal_charge,
            chiral_tag=chirality_encoding.get(chiral_tag_encoded, rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED),
            hybridization=hybridization_encoding[ hybridization_encoded],
            num_explicit_hs=num_explicit_hs,
            is_aromatic=is_aromatic,
            isotope=isotope,
            AtomName = AtomName,
            ResID=ResID,

        )



    for (bond_type_encoded, is_peptide_bond, bond_start, bond_end) in mol_json['edges_tuple']:
        G.add_edge(
            bond_start,
            bond_end,
            bond_type=bond_type_encoding[bond_type_encoded],
            is_peptide_bond = is_peptide_bond,)
    return G

def json_to_mol(mol_json:dict):
    G_mol = json_to_nx(mol_json)
    return nx_to_mol(G_mol)


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