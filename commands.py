"""
#**pepseq.commands**

Provide several sample math calculations.

This module allows the user to make mathematical calculations.


Examples:
    >>> from pepseq.commands import pepseq_to_smiles, calculate_json_from,
        read_smiles, augment_db_json_command
    >>> pepseq_to_smiles('CH3~SCAFC~NH2')
    
    >>> calculate_json_from('CH3~SC{R1}AFC~NH2', '[*1]CCC') calculations.multiply(2.0, 4.0)
    
    >>> read_smiles('mypeptide.smi', 'myppeptide_out')

    >>> augment_db_json_command('my_monomers.sdf', 'augmented_db.json')

The module contains the following functions:

- `pepseq_to_smiles(pepseq, out, db_path)` - Returns the SMILES code for pepseq.
- `calculate_json_from(sequence, mod_smiles, out, db_path)` - Returns JSON dict 
   containing info about amino acid sequence and sequence modifications
- `read_smiles(smiles_filename, out, db_path, v) - Reads SMILES from file  and
   Writes parsed Pepseq and its modifications into separate txt files.
- `augment_db_json_command(sdf_path, out)` - Reads additional (Modified) Peptide
   building blocks from SDF file. Adds them to database and outputs it to file.

"""

import json
import pkgutil
from tqdm import tqdm

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
        db_path: str = db_path,
        v: bool = False
        ):
    """
    Input SMILES filepath 


    python3.10 commands.py read-smiles 'stuff_in.smi'   --db-path augmented_db.json --out stuff

    python3.10 commands.py read-smiles 'smilesy.smi' --out smilesy_out
    
    ls 
    """
    with open(smiles_filename) as fp:
        s = fp.read()
    smiles_list = [smiles for smiles in s.split('\n') if smiles]

    pepseq_list = []

    for smiles in tqdm(smiles_list):
        kwargs = {}

        if db_path is not None:
            with open(db_path) as fp:
                db_json = json.load(fp)
                kwargs['db_json'] = db_json


        pepseq_format, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
            smiles, **kwargs
        )
        if v:
            print('Sequence in pepseq format: ', pepseq_format)
            print('List of modification by SMILES codes with points of attachment: ', mod_smiles)
        pepseq_list.append( pepseq_format )
        if type(mod_smiles) == list:
            if mod_smiles:
                mod_smiles_list.append( '\t'.join(mod_smiles) )
            else:
                mod_smiles_list.append('NO (NON TERMINAL) SEQUENCE MODIIFICATIONS PRESENT')
        else:
            mod_smiles_list.append( mod_smiles )
    
    out_pepseq = '%s.pepseq' %(out)
    out_mod_smiles = '%s.smi' %(out)

    with open(out_pepseq, 'w') as fp:
        fp.write( '\n'.join( pepseq_list) )
        fp.flush()
    
    if type(mod_smiles) == str:
        mod_smiles_list = [mod_smiles]

    mod_smiles_list_str = '\n'.join(mod_smiles_list)

    with open(out_mod_smiles, 'w') as fp:
        fp.write(mod_smiles_list_str)
        fp.flush()

    return


@app.command()
def augment_db_json_command(sdf_path='sdf_file.sdf',
                            out='db_json_augmented.json'):
    """

    Input:
        sdf_path: path for SDF file containing info on (Modified) Peptide 
        building blocks
    
    Output:
        out: path for db.json dictionary containing databse unriched with
          (Modified) Peptide building blocks declared in SDF file
    
    Action:
        SDF file is read by rdkit. Molecules contained in SDF file are processed
        default exit atom (atom to attach modifications unless
        specified otherwise) is computed by set_default_exit_atom_function
        from this CXSMILES and CXSMARTS are created together with SMILES code
        depicting basic molecule (without radicals attached).

        Modified Database JSON is written to output file

    """

    df_sdf = rdkit.Chem.PandasTools.LoadSDF(sdf_path)

    db_json_augmented = augment_db_json(db_json, df_sdf=df_sdf, name_column = 'm_abbr', mol_colname='ROMol')

    with open(out, 'w') as fp:
        json.dump(db_json_augmented, fp)

    return


if __name__ == "__main__":
    app()
