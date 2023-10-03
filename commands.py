import json
import pkgutil

import typer
from pepseq.read import from_pepseq
from pepseq.functions import calculate

from typing import List, Optional
from typing_extensions import Annotated

import rdkit.Chem.PandasTools

from pepseq.BuildPeptideJSONFromSMILES import \
    from_smiles_to_pepseq_and_one_mod_smiles_strings
from augmenting_db_json import augment_db_json


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


app = typer.Typer()


@app.command()
def pepseq_to_smiles(
    pepseq: str,
    out: str = None,
    db_path: str = None
) -> str:
    """
    Input Sequence in Pepseq format and generate molecule SMILES code.
    Example Pepseq Format Sequence:

    "[CSPS[1*]]~Y{aMeAla}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[1*]]"

     "[CSPS[1*]]~Y{aMeAla}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[1*]]"

    python3.10 commands.py pepseq-to-smiles  "[CSPS[1*]]~Y{aMeAla}QGT{Orn}FTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[1*]]"  --db-path augmented_db.json --out out.smi
    """
    if db_path is not None:
        with open(db_path) as fp:
            db_json = json.load(fp)

    peptide = from_pepseq(pepseq, db_json=db_json)
    if out is not None:
        with open(out,'w') as fp:
            fp.write(peptide.smiles)
            fp.flush()

    return peptide.smiles


@app.command()
def calculate_json_from(
    sequence: str,
    mod_smiles: Annotated[Optional[List[str]], typer.Option(default=None)] = None,
    out: str = None,
    db_path: str = None
    ):
    """

    Calculate peptide json from sequence and list of modification SMILES

    Example:

    python3.10 commands.py calculate-json-from
       '{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF'
         --mod-smiles '[1*]CNCC[2*]' --mod-smiles '[3*]CNCCSP'

    """
    kwargs = {}

    if db_path is not None:
        with open(db_path) as fp:
            db_json = json.load(fp)
            kwargs['db_json'] = db_json

    args = (sequence, mod_smiles)
    result = calculate(*args, **kwargs)
    print('complete_smiles: ', result.get('complete_smiles'))
    print('sequence: ', result.get('sequence'))
    if out is not None:
        with open(out,'w') as fp:
            json.dump(result, fp)
    return result



@app.command()
def read_smiles(
        smiles_filename: str,
        out: str = 'out',
        db_path: str = None
        ):
    """
    python3.10 commands.py read-smiles --smiles-filename  augmented.smi  --db-path augmented_db.json --out stuff

    """
    with open(smiles_filename) as fp:
        smiles = fp.readline().strip()

    kwargs = {}

    if db_path is not None:
        with open(db_path) as fp:
            db_json = json.load(fp)
            kwargs['db_json'] = db_json


    pepseq_format, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        smiles, **kwargs
    )    

    print('Sequence in pepseq format: ', pepseq_format)
    print('List of modification by SMILES codes with points of attachment: ', mod_smiles)
    

    out_pepseq = '%s.pepseq' %(out)
    out_mod_smiles = '%s.smi' %(out)

    with open(out_pepseq, 'w') as fp:
        fp.write(pepseq_format)
        fp.flush()
    
    if type(mod_smiles) == str:
        mod_smiles = [mod_smiles]

    mod_smiles_list_str = '\n'.join(mod_smiles)

    with open(out_mod_smiles, 'w') as fp:
        fp.write(mod_smiles_list_str)
        fp.flush()

    return


@app.command()
def augment_db_json_command(sdf_path='sdf_file.sdf',
                            out='db_json_augmented.json'):
    df_sdf = rdkit.Chem.PandasTools.LoadSDF(sdf_path)

    db_json_augmented = augment_db_json(db_json, df_sdf=df_sdf, name_column = 'm_abbr', mol_colname='ROMol')

    with open(out, 'w') as fp:
        json.dump(db_json_augmented, fp)

    return


if __name__ == "__main__":
    app()
