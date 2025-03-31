import json
import copy
import pandas as pd
import rdkit
import rdkit.Chem
from pepseq.augmenting_db_json import augment_db_json, get_smiles_set
from pepseq.tests_new.helpers import check_db_json_all_keys


db_json_ori_CDG = {
  'coding': {
    'l_proteogenic_3letter': {'Cys': 'C', 'Asp': 'D', 'Gly': 'G'},
    'd_proteogenic_3letter': {'cys': 'C', 'asp': 'D', 'gly': 'G'},
    'aa_codes': {'C': 'CYS', 'D': 'ASP', 'G': 'GLY'},
    'd_proteogenic_2letter': {'dC': 'c', 'dD': 'd', 'dG': 'g'},
    'd_proteogenic_4letter': {'dCys': 'c', 'dAsp': 'd', 'dGly': 'g'},
    'modified_aa_codes': {
      'aMeAla': 'X',
      '2-Aminoisobutyric acid': 'X',
      '2-Methylalanine': 'X'
      },
    'modified_aa_codes_reverse': {'X': 'aMeAla'}
    },
  'physico_chemical_properties': {
    'hydrophobicity_index_pH2': {'C': 53, 'G': 0, 'D': -18},
    'hydrophobicity_index_pH7': {'C': 49, 'G': 0, 'D': -55},
    'acidic_aa': {'D': ''},
    'basic_aa': {},
    'protein_weights': {'C': 121.1582, 'D': 133.1027, 'G': 75.0666}
    },
  'l_proteogenic': ['C', 'D', 'G'],
  'd_proteogenic': ['c', 'd', 'g',],
  'aa_order': ['g', 'G', 'dD', 'C_nonmod', 'C', 'd', 'D'],
  'n_terms_order': ['H', 'Ac', 'Pro', 'Myr', 'Palm'],
  'c_terms_order': ['OH', 'NH2', 'N-Me'],
  'protein_letters': 'CDG',
  'smiles': {
    'n_terms': {
      'CH3': {
        'smiles': 'C',
        'smiles_radical': '[*:1]C',
        'three_letter_code': 'CH3',
        'smarts': 'C',
        'smarts_comments': ''
        },
      },
    'c_terms': {
      'NH2': {
        'smiles': '[NH2]',
        'smiles_radical': '[*:1]N',
        'three_letter_code': 'NH2',
        'smarts': '[NH2]',
        'smarts_comments': ''
        },
      },
    'aa': {
      'C': {
        'smiles': 'N[C@H](C=O)CS',
        'smiles_radical': '[*:1]N[C@@H](CS)C([*:2])=O |atomProp:0.dummyLabel.1*:4.AtomName.SG:4.atomLabel.SG:6.dummyLabel.2*|',
        'three_letter_code': 'Cys',
        'smarts': '[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)][C&H1&X4]([C&H2&X4][S&X2&H1,S&X1&H0&-,S&X2])[C](=[O&X1]) |atomProp:3.AtomName.SG|',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'C_nonmod': {
        'smiles': 'N[C@H](C=O)CS',
        'smiles_radical': 'N([*:1])[C@H](C([*:2])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|',
        'three_letter_code': 'Cys',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'D': {
        'smiles': 'N[C@H](C=O)CC(=O)O',
        'smiles_radical': 'N([*:1])[C@H](C([*:2])=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|',
        'three_letter_code': 'Asp',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])',
        'smarts_comments': 'Hits acid and conjugate base. Also hits Glu side chain when used alone.'
        },
      'G': {
        'smiles': 'NCC=O',
        'smiles_radical': '[*:1]NCC([*:2])=O |$;;CA;;;$,atomProp:0.dummyLabel.1*:2.AtomName.CA:4.dummyLabel.2*|',
        'three_letter_code': 'Gly',
        'smarts': 'N[CX4H2][CX3](=[OX1])',
        'smarts_comments': ''
        },
      'g': {
        'smiles': 'NCC=O',
        'smiles_radical': 'N([*:1])CC([*:2])=O |$_N;_R1;_CA;_CO;_R2$|',
        'three_letter_code': 'Gly',
        'smarts': 'N[CX4H2][CX3](=[OX1])[O,N]',
        'smarts_comments': ''
        },
      'c': {
        'smiles': 'N[C@@H](C=O)CS',
        'smiles_radical': 'N([*:1])[C@@H](C([*:2])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|',
        'three_letter_code': 'cys',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'c_nonmod': {
        'smiles': 'N[C@@H](C=O)CS',
        'smiles_radical': 'N([*:1])[C@@H](C([*:2])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|',
        'three_letter_code': 'cys',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'd': {
        'smiles': 'N[C@@H](C=O)CC(=O)O',
        'smiles_radical': 'N([*:1])[C@@H](C([*:2])=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|',
        'three_letter_code': 'asp',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])',
        'smarts_comments': 'Hits acid and conjugate base. Also hits Glu side chain when used alone.'
        },
      'dD': {
        'smiles': 'N[C@@H](C=O)CC=O',
        'smiles_radical': '[*:1]N[C@H](CC=O)C([*:2])=O |$_AV:*;;;;;;;*;$,atomProp:0.hybridization.UNSPECIFIED:0.dummyLabel.R1:0.num_explicit_hs.0:1.hybridization.SP3:1.num_explicit_hs.0:2.hybridization.SP3:2.num_explicit_hs.1:3.hybridization.SP3:3.num_explicit_hs.0:4.hybridization.SP2:4.num_explicit_hs.0:4.DefaultExitAtom.True:5.hybridization.SP2:5.num_explicit_hs.0:6.hybridization.SP2:6.num_explicit_hs.0:7.hybridization.UNSPECIFIED:7.dummyLabel.R2:7.num_explicit_hs.0:8.hybridization.SP2:8.num_explicit_hs.0|',
        'three_letter_code': 'dD',
        'smarts': '[#7]-[#6@&H1](-[#6]-[#6]=[#8])-[#6]=[#8] |atomProp:0.hybridization.SP3:0.num_explicit_hs.0:1.hybridization.SP3:1.num_explicit_hs.1:2.hybridization.SP3:2.num_explicit_hs.0:3.hybridization.SP2:3.num_explicit_hs.0:4.hybridization.SP2:4.num_explicit_hs.0:5.hybridization.SP2:5.num_explicit_hs.0:6.hybridization.SP2:6.num_explicit_hs.0|',
        'smarts_comments': ''
        },
      }
    }
  }


db_json_CDG_augmented_fixture = {
  'coding': {
    'l_proteogenic_3letter': {'Cys': 'C', 'Asp': 'D', 'Gly': 'G'},
    'd_proteogenic_3letter': {'cys': 'C', 'asp': 'D', 'gly': 'G'},
    'aa_codes': {'C': 'CYS', 'D': 'ASP', 'G': 'GLY'},
    'd_proteogenic_2letter': {'dC': 'c', 'dD': 'd', 'dG': 'g'},
    'd_proteogenic_4letter': {'dCys': 'c', 'dAsp': 'd', 'dGly': 'g'},
    'modified_aa_codes': {
      'aMeAla': 'X',
      '2-Aminoisobutyric acid': 'X',
      '2-Methylalanine': 'X'
      },
    'modified_aa_codes_reverse': {'X': 'aMeAla'}
    },
  'physico_chemical_properties': {
    'hydrophobicity_index_pH2': {'C': 53, 'G': 0, 'D': -18},
    'hydrophobicity_index_pH7': {'C': 49, 'G': 0, 'D': -55},
    'acidic_aa': {'D': ''},
    'basic_aa': {},
    'protein_weights': {'C': 121.1582, 'D': 133.1027, 'G': 75.0666}
    },
  'l_proteogenic': ['C', 'D', 'G'],
  'd_proteogenic': ['c', 'd', 'g'],
  'aa_order': ['g', 'G', 'dD', 'C_nonmod', 'C', 'd', 'D'],
  'n_terms_order': ['H', 'Ac', 'Pro', 'Myr', 'Palm'],
  'c_terms_order': ['OH', 'NH2', 'N-Me'],
  'protein_letters': 'CDG',
  'smiles': {
    'n_terms': {
      'CH3': {
        'smiles': 'C',
        'smiles_radical': '[*:1]C',
        'three_letter_code': 'CH3',
        'smarts': 'C',
        'smarts_comments': ''
        }
      },
    'c_terms': {
      'NH2': {
        'smiles': '[NH2]',
        'smiles_radical': '[*:1]N',
        'three_letter_code': 'NH2',
        'smarts': '[NH2]',
        'smarts_comments': ''
        }
      },
    'aa': {
      'C': {
        'smiles': 'N[C@H](C=O)CS',
        'smiles_radical': '[*:1]N[C@@H](CS)C([*:2])=O |atomProp:0.dummyLabel.1*:4.AtomName.SG:4.atomLabel.SG:6.dummyLabel.2*|',
        'three_letter_code': 'Cys',
        'smarts': '[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)][C&H1&X4]([C&H2&X4][S&X2&H1,S&X1&H0&-,S&X2])[C](=[O&X1]) |atomProp:3.AtomName.SG|',
        'smarts_comments': 'Hits acid and conjugate base'},
      'C_nonmod': {
        'smiles': 'N[C@H](C=O)CS',
        'smiles_radical': 'N([*:1])[C@H](C([*:2])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|',
        'three_letter_code': 'Cys',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'D': {
        'smiles': 'N[C@H](C=O)CC(=O)O',
        'smiles_radical': 'N([*:1])[C@H](C([*:2])=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|',
        'three_letter_code': 'Asp',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])',
        'smarts_comments': 'Hits acid and conjugate base. Also hits Glu side chain when used alone.'
        },
      'G': {
        'smiles': 'NCC=O',
        'smiles_radical': '[*:1]NCC([*:2])=O |$;;CA;;;$,atomProp:0.dummyLabel.1*:2.AtomName.CA:4.dummyLabel.2*|',
        'three_letter_code': 'Gly',
        'smarts': 'N[CX4H2][CX3](=[OX1])',
        'smarts_comments': ''
        },
      'g': {
        'smiles': 'NCC=O',
        'smiles_radical': 'N([*:1])CC([*:2])=O |$_N;_R1;_CA;_CO;_R2$|',
        'three_letter_code': 'Gly',
        'smarts': 'N[CX4H2][CX3](=[OX1])[O,N]',
        'smarts_comments': ''
        },
      'c': {
        'smiles': 'N[C@@H](C=O)CS',
        'smiles_radical': 'N([*:1])[C@@H](C([*:2])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|',
        'three_letter_code': 'cys',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'c_nonmod': {
        'smiles': 'N[C@@H](C=O)CS',
        'smiles_radical': 'N([*:1])[C@@H](C([*:2])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|',
        'three_letter_code': 'cys',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]',
        'smarts_comments': 'Hits acid and conjugate base'
        },
      'd': {
        'smiles': 'N[C@@H](C=O)CC(=O)O',
        'smiles_radical': 'N([*:1])[C@@H](C([*:2])=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|',
        'three_letter_code': 'asp',
        'smarts': '[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])',
        'smarts_comments': 'Hits acid and conjugate base. Also hits Glu side chain when used alone.'
        },
      'dD': {
        'smiles': 'N[C@@H](C=O)CC=O',
        'smiles_radical': '[*:1]N[C@H](CC=O)C([*:2])=O |$_AV:*;;;;;;;*;$,atomProp:0.hybridization.UNSPECIFIED:0.dummyLabel.R1:0.num_explicit_hs.0:1.hybridization.SP3:1.num_explicit_hs.0:2.hybridization.SP3:2.num_explicit_hs.1:3.hybridization.SP3:3.num_explicit_hs.0:4.hybridization.SP2:4.num_explicit_hs.0:4.DefaultExitAtom.True:5.hybridization.SP2:5.num_explicit_hs.0:6.hybridization.SP2:6.num_explicit_hs.0:7.hybridization.UNSPECIFIED:7.dummyLabel.R2:7.num_explicit_hs.0:8.hybridization.SP2:8.num_explicit_hs.0|',
        'three_letter_code': 'dD',
        'smarts': '[#7]-[#6@&H1](-[#6]-[#6]=[#8])-[#6]=[#8] |atomProp:0.hybridization.SP3:0.num_explicit_hs.0:1.hybridization.SP3:1.num_explicit_hs.1:2.hybridization.SP3:2.num_explicit_hs.0:3.hybridization.SP2:3.num_explicit_hs.0:4.hybridization.SP2:4.num_explicit_hs.0:5.hybridization.SP2:5.num_explicit_hs.0:6.hybridization.SP2:6.num_explicit_hs.0|',
        'smarts_comments': ''
        },
      'A': {
        'smiles': 'C[C@H](N)C=O',
        'smiles_radical': '*N[C@@H](C)C(*)=O |$_AV:*;;;;;*;$,atomProp:0.dummyLabel.R1:0.hybridization.UNSPECIFIED:0.num_explicit_hs.0:1.num_explicit_hs.0:1.hybridization.SP3:2.num_explicit_hs.1:2.hybridization.SP3:3.num_explicit_hs.0:3.hybridization.SP3:4.num_explicit_hs.0:4.hybridization.SP2:5.dummyLabel.R2:5.hybridization.UNSPECIFIED:5.num_explicit_hs.0:6.num_explicit_hs.0:6.hybridization.SP2|',
        'three_letter_code': 'A',
        'smarts': '[#0]-[#7]-[#6@@H](-[#6])-[#6](-[#0])=[#8] |$_AV:*;;;;;*;$,atomProp:0.dummyLabel.R1:0.hybridization.UNSPECIFIED:0.num_explicit_hs.0:1.num_explicit_hs.0:1.hybridization.SP3:2.num_explicit_hs.1:2.hybridization.SP3:3.num_explicit_hs.0:3.hybridization.SP3:4.num_explicit_hs.0:4.hybridization.SP2:5.dummyLabel.R2:5.hybridization.UNSPECIFIED:5.num_explicit_hs.0:6.num_explicit_hs.0:6.hybridization.SP2|',
        'smarts_comments': ''
        },
      'L': {
        'smiles': 'CC(C)C[C@H](N)C=O',
        'smiles_radical': '*N[C@@H](CC(C)C)C(*)=O |$_AV:*;;;;;;;;*;$,atomProp:0.dummyLabel.R1:0.hybridization.UNSPECIFIED:0.num_explicit_hs.0:1.num_explicit_hs.0:1.hybridization.SP3:2.num_explicit_hs.1:2.hybridization.SP3:3.num_explicit_hs.0:3.hybridization.SP3:4.num_explicit_hs.0:4.hybridization.SP3:5.num_explicit_hs.0:5.hybridization.SP3:6.num_explicit_hs.0:6.hybridization.SP3:7.num_explicit_hs.0:7.hybridization.SP2:8.dummyLabel.R2:8.hybridization.UNSPECIFIED:8.num_explicit_hs.0:9.num_explicit_hs.0:9.hybridization.SP2|',
        'three_letter_code': 'L',
        'smarts': '[#0]-[#7]-[#6@@H](-[#6]-[#6](-[#6])-[#6])-[#6](-[#0])=[#8] |$_AV:*;;;;;;;;*;$,atomProp:0.dummyLabel.R1:0.hybridization.UNSPECIFIED:0.num_explicit_hs.0:1.num_explicit_hs.0:1.hybridization.SP3:2.num_explicit_hs.1:2.hybridization.SP3:3.num_explicit_hs.0:3.hybridization.SP3:4.num_explicit_hs.0:4.hybridization.SP3:5.num_explicit_hs.0:5.hybridization.SP3:6.num_explicit_hs.0:6.hybridization.SP3:7.num_explicit_hs.0:7.hybridization.SP2:8.dummyLabel.R2:8.hybridization.UNSPECIFIED:8.num_explicit_hs.0:9.num_explicit_hs.0:9.hybridization.SP2|',
        'smarts_comments': ''
        }
      }
    }
  }

df_sdf_mask_copy_read = pd.read_csv('data/df_sdf_mask_copy.csv')
df_sdf_mask_copy_read['ROMol'] = df_sdf_mask_copy_read['ROMol_CXSmiles'].apply(
    lambda x: rdkit.Chem.MolFromSmiles(x))


smiles_set_fixture = set(
  [
    'NCC=O',
    'N[C@@H](C=O)CC(=O)O',
    'N[C@@H](C=O)CC=O',
    'N[C@@H](C=O)CS',
    'N[C@H](C=O)CC(=O)O',
    'N[C@H](C=O)CS'
    ]
  )


smiles_set_fixture2 = set(
  [
    'NCC=O',
    'N[C@@H](C=O)CC(=O)O',
    'N[C@@H](C=O)CC=O',
    'N[C@@H](C=O)CS',
    'N[C@H](C=O)CC(=O)O',
    'N[C@H](C=O)CS'
    ]
  )


def test_augment_db_json():
    """
    Tests the augment_db_json function by performing the following steps:
    1. Creates a deep copy of the original database JSON (db_json_ori_CDG).
    2. Augments the copied JSON using the augment_db_json function with specified parameters.
    3. Writes the augmented JSON to a file named 'db_json_CDG_augmented.json'.
    4. Asserts that the augmented JSON matches the expected fixture (db_json_CDG_augmented_fixture).
    This test ensures that the augment_db_json function correctly processes and augments the input JSON data.
    """
    db_json_ori_CDG_copy = copy.deepcopy(db_json_ori_CDG)
    db_json_CDG_augmented = augment_db_json(
        db_json_ori_CDG_copy, df_sdf=df_sdf_mask_copy_read, name_column="m_abbr", mol_colname="ROMol"
    )
    assert check_db_json_all_keys(db_json_CDG_augmented, db_json_CDG_augmented_fixture)
    with open('data/db_json_CDG_augmented.json', 'w') as fp:
        json.dump(db_json_CDG_augmented, fp)

    return


def test_get_smiles_set():
    """
    Test the get_smiles_set function to ensure it returns the correct set of SMILES strings.
    This test creates a deep copy of the original database JSON object (db_json_ori_CDG),
    passes it to the get_smiles_set function, and asserts that the returned set of SMILES
    strings matches the expected fixture (smiles_set_fixture2).
    Returns:
      None
    """
    db_json_ori_CDG_copy = copy.deepcopy(db_json_ori_CDG)
    smiles_set = get_smiles_set(db_json_ori_CDG_copy)
    assert smiles_set == smiles_set_fixture2
    return
