from commands import read_smiles, get_ketcher_param, calculate_json_from


def test_read_smiles():
    #pepseq, smiles = data_pepseq_smiles
    peptide = read_smiles('outpath.smi')
    #assert smiles_are_identical(peptide.smiles, smiles)


# def get_ketcher_param(smi):

def test_get_ketcher_param():

    return


def test_calculate_json_from():
    #pepseq, smiles = data_pepseq_smiles
    calculate_json_from('CH3~SC{Cys(R1)}AFC~NH2', '[*:1]CCC', ketcher=True)
    #assert smiles_are_identical(peptide.smiles, smiles)
