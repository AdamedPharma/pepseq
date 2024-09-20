import json
import pkgutil

from pepseq.BuildPeptideJSONFromSMILES import  from_smiles_to_pepseq_and_mod_smiles_strings, \
    decompose_peptide_smiles_with_termini, decompose_peptide_smiles, from_smiles_to_pepseq_and_mod_smiles_strings, \
    from_smiles_to_pepseq_and_one_mod_smiles_strings



#from_smiles_to_pepseq_and_mod_smiles_strings(smiles: str, db_json: dict, n_subst_limit=None)

#    peptide_json = decompose_peptide_smiles_with_termini(smiles, db_json, n_subst_limit=n_subst_limit)
#    pepseq_format = peptide_json["pepseq_format"]
#    mod_smiles_list = [
#        ext_mod["smiles"] for ext_mod in peptide_json["external_modifications"]
#    ]

db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)



def test_from_smiles_to_pepseq_and_mod_smiles_strings():
    smiles = '[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCCSP)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O'
    n_subst_limit = None

    pepseq, mod_smiles = from_smiles_to_pepseq_and_mod_smiles_strings(
        smiles, db_json, n_subst_limit=n_subst_limit, ketcher=True
        )
    
    assert pepseq == 'H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH'
    assert mod_smiles == ['C(C[*:2])NC[*:1]', 'PSCCNC[*:3]']

    pepseq, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(smiles, db_json, n_subst_limit=None, ketcher=True)

    assert pepseq == 'H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH'
    assert mod_smiles == ['C(C[*:2])NC[*:1]', 'PSCCNC[*:3]']




def test_decompose_peptide_smiles():
    smiles_fixture = '[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCCSP)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O'
    peptide_json_fixture = {
        'sequence': 'CACDAPEPsEQCGCDEF',
        'internal_modifications': [],
        'external_modifications': [
            {
                'smiles': 'C(C[*:2])NC[*:1]',
                'max_attachment_point_id': 2,
                'attachment_points_on_sequence': {
                    1: {
                        'attachment_point_id': 1,
                        'ResID': '1',
                        'AtomName': 'SG',
                        'ResidueName': ''
                        },
                    2: {
                        'attachment_point_id': 2,
                        'ResID': '12',
                        'AtomName': 'SG',
                        'ResidueName': ''
                        }
                    }
                },
            {
                'smiles': 'PSCCNC[*:1]',
                'max_attachment_point_id': 1,
                'attachment_points_on_sequence': {
                    1: {
                        'attachment_point_id': 1,
                        'ResID': '14',
                        'AtomName': 'SG',
                        'ResidueName': ''
                        }
                    }
                },
            {
                'smiles': 'O[*:1]',
                'max_attachment_point_id': 1,
                'attachment_points_on_sequence': {
                    1: {
                        'attachment_point_id': 1,
                        'ResID': '17',
                        'AtomName': 'CO',
                        'ResidueName': ''
                        }
                    }
                }
            ]
        }
    
    
    peptide_json_w_termini_fixture = {'sequence': 'CACDAPEPsEQCGCDEF',
 'internal_modifications': [],
 'external_modifications': [{'smiles': 'C(C[*:2])NC[*:1]',
   'max_attachment_point_id': 2,
   'attachment_points_on_sequence': {2: {'attachment_point_id': 2,
     'ResID': '12',
     'AtomName': 'SG',
     'ResidueName': ''},
    1: {'attachment_point_id': 1,
     'ResID': '1',
     'AtomName': 'SG',
     'ResidueName': ''}}},
  {'smiles': 'PSCCNC[*:3]',
   'max_attachment_point_id': 3,
   'attachment_points_on_sequence': {3: {'attachment_point_id': 1,
     'ResID': '14',
     'AtomName': 'SG',
     'ResidueName': ''}}}],
 'C_terminus': 'OH',
 'N_terminus': 'H',
 'pepseq_format': 'H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH'}
    n_subst_limit = None
    peptide_json = decompose_peptide_smiles(smiles_fixture, db_json, n_subst_limit=n_subst_limit, ketcher = True)

    assert peptide_json == peptide_json_fixture
    peptide_json_w_termini = decompose_peptide_smiles_with_termini(smiles_fixture, db_json, n_subst_limit=n_subst_limit, ketcher = True)
    assert peptide_json_w_termini == peptide_json_w_termini_fixture

    return