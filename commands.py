"""
#**pepseq.commands**

Provide several sample math calculations.

Pepseq: chemoinformatic library to work flawlessly with modified peptide structures.
Pepseq is a library to manage and analyse modified peptide sequences.
The pepseq format bridges the gap between simplicity of FASTA and expressiveness
 of SMILES for modified peptides, addressing a real need in chemoinformatics. 

This module allows user to enter modified peptide structure in SMILES format,
output it in pepseq_format.

Pepseq format is composed of one (or two) parts. First part is a sequence in FASTA
 format with annotated N and C termini. The second (optional) part is modification(s)
 defined in SMILES format (with R group(s) in atoms connecting sequence
 to modification(s)). Connecting residues on peptide sequence are also annotated.  
Reverse transformation is also possible.
Monomers used in sequence are defined in database in JSON format (db.json).
It is possible to augment database with new monomers, e.g. defined in SDF file. 

Examples:
    >>> from pepseq.commands import pepseq_to_smiles, calculate_json_from,
        read_smiles, augment_db_json_command
    >>> pepseq_to_smiles('CH3~SCAFC~NH2')

    >>> calculate_json_from('CH3~SC{R1}AFC~NH2', '[*:1]CCC')

    >>> read_smiles('mypeptide.smi', 'myppeptide_out')

    >>> augment_db_json_command('my_monomers.sdf', 'augmented_db.json')

The module contains the following functions:

- `pepseq_to_smiles(pepseq, out, db_path)` - Returns the SMILES code for pepseq.
- `calculate_json_from(sequence, mod_smiles, out, db_path)` - Returns JSON dict
   containing info about amino acid sequence and sequence modifications
- `read_smiles(smiles_filename, out, db_path, v)` - Reads SMILES from file and
   Writes parsed Pepseq and its modifications into separate txt files.
- `augment_db_json_command(sdf_path, out)` - Reads additional (Modified) Peptide
   building blocks from SDF file. Adds them to database and outputs it to file.

"""

import json
import os

from tqdm import tqdm

import typer
from pepseq.read import from_pepseq
from pepseq.functions import calculate

from typing import List, Optional, Union
from typing_extensions import Annotated

import rdkit.Chem.PandasTools

from pepseq.BuildPeptideJSONFromSMILES import (
    from_smiles_to_pepseq_and_one_mod_smiles_strings,
)
from pepseq.augmenting_db_json import augment_db_json
from pepseq.Peptide.database.db_functions import load_db_json

dir_path = os.path.dirname(__file__)
db_path = os.path.join(dir_path, "pepseq/Peptide/database/db.json")
with open(db_path) as fp:
    db_json = json.load(fp)


app = typer.Typer()


@app.command()
def pepseq_to_smiles(pepseq: str, out: Union[str, None] = None, db_path: Union[str, None] = None) -> str:
    """
    Return SMILES code for Modified Peptide given by text in Pepseq Format.

    :param pepseq: string in pepseq format.
    :type pepseq: str or None
    :param out: output path.
    :type out: str or None
    :param db_path: Optional db path.
    :type db_path: str or None

    :return: SMILES code.
    :rtype: str

    Example Pepseq Format Sequence:

    "[CSPS[*:1]]~Y{aMeAla}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[*:1]]"

     "[CSPS[*:1]]~Y{aMeAla}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[*:1]]"

    python commands.py pepseq-to-smiles
        "[CSPS[*:1]]~Y{aMeAla}QGT{Orn}FTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[*:1]]"
            --db-path augmented_db.json --out out.smi
    """
    if db_path is None:
        db_path = os.path.join(dir_path, "pepseq/Peptide/database/db.json")
    with open(db_path) as fp:
        db_json = json.load(fp)

    peptide = from_pepseq(pepseq, db_json=db_json)
    if out is not None:
        with open(out, "w") as fp:
            fp.write(peptide.smiles)
            fp.flush()

    return peptide.smiles


@app.command()
def calculate_json_from(
    sequence: str,
    mod_smiles: Annotated[
        Optional[Union[List[str],None]], "List of modification SMILES codes"
    ] = None,
    out: str = None,
    db_path: str = None,
) -> dict:
    """

    Calculate peptide json from sequence and list of modification SMILES

    :param sequence: string in pepseq format.
    :type sequence: str or None
    :param mod_smiles: string list of SMILES codes.
    :type mod_smiles: list[str] or None

    :param out: output path.
    :type out: str or None
    :param db_path: Optional db path.
    :type db_path: str or None

    :return: SMILES code.
    :rtype: dict

    Example:

    import pepseq
    pepseq.commands.calculate_json_from
    '{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF'
    --mod-smiles '[*:1]CNCC[*:2]' --mod-smiles '[*:3]CNCCSP'

    """
    kwargs = {}

    db_json = load_db_json(db_path)
    if db_json is not None:
        kwargs["db_json"] = db_json

    if type(mod_smiles) == str:
        mod_smiles = [mod_smiles]

    args = (sequence, mod_smiles)
    result = calculate(*args, **kwargs)
    print("complete_smiles: ", result.get("complete_smiles"))
    print("sequence: ", result.get("sequence"))
    if out is not None:
        with open(out, "w") as fp:
            json.dump(result, fp)
    return result


@app.command()
def read_smiles(
    smiles_filename: str,
    out: str = "out",
    db_path: str = db_path,
    v: bool = False,
) -> list:
    """

    Calculate sequence(s) in Pepseq Format and  list of modification SMILES from
    Modified Peptide structure given in SMILES.
    Output results into files

    :param smiles_filename: string in pepseq format.
    :type smiles_filename: str or None
    :param mod_smiles: string list of SMILES codes.
    :type mod_smiles: list[str] or None

    :param out: output path.
    :type out: str or None
    :param db_path: Optional db path.
    :type db_path: str or None
    :param v: bool switch regulating verbosity.
    :type v: bool or None

    :return: result_list: list of results given as list of sequences in Pepseq Format. and list of modification SMILES
    :rtype: list

    Example:

    python commands.py read-smiles 'stuff_in.smi'   --db-path augmented_db.json --out stuff

    python commands.py read-smiles 'smilesy.smi' --out smilesy_out

    """
    with open(smiles_filename) as fp:
        s = fp.read()

    smiles_list = [smiles for smiles in s.split("\n") if smiles]

    pepseq_list = []

    for smiles in tqdm(smiles_list):
        kwargs = {}
        db_json  = load_db_json(db_path)
        if db_json is not None:
            kwargs["db_json"] = db_json

        mod_smiles_list = []

        pepseq_format, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
            smiles, **kwargs
        )

        if v:
            print("Sequence in pepseq format: ", pepseq_format)
            print(
                "List of modification by SMILES codes with points of attachment: ",
                mod_smiles,
            )
        pepseq_list.append(pepseq_format)
        if type(mod_smiles) == list:
            if mod_smiles:
                mod_smiles_list.append("\t".join(mod_smiles))
            else:
                mod_smiles_list.append(
                    "NO (NON TERMINAL) SEQUENCE MODIIFICATIONS PRESENT"
                )
        else:
            mod_smiles_list.append(mod_smiles)

    out_pepseq = f"{out}.pepseq"

    out_mod_smiles = f"{out}.smi"

    with open(out_pepseq, "w") as fp:
        fp.write("\n".join(pepseq_list))
        fp.flush()

    if type(mod_smiles) == str:
        mod_smiles_list = [mod_smiles]

    mod_smiles_list_str = "\n".join(mod_smiles_list)

    with open(out_mod_smiles, "w") as fp:
        fp.write(mod_smiles_list_str)
        fp.flush()

    return [pepseq_list, mod_smiles_list]


@app.command()
def augment_db_json_command(
    sdf_path="sdf_file.sdf", out="db_json_augmented.json"
) -> dict:
    """
    SDF file is read by rdkit. Molecules contained in SDF file are processed
        default exit atom (atom to attach modifications unless
        specified otherwise) is computed by set_default_exit_atom_function
        from this CXSMILES and CXSMARTS are created together with SMILES code
        depicting basic molecule (without radicals attached).

    :param sdf_path: path for SDF format file containing new monomers to be
        inserted into database.
    :type sdf_path: str or None

    :param out: output path for db_json copy with monomers and modifications
        added from SFD file.
    :type out: str or None

    :return: db_json database copy with augmented by monomers from SDF file.
    :rtype: dict
    """

    df_sdf = rdkit.Chem.PandasTools.LoadSDF(sdf_path)

    db_json_augmented = augment_db_json(
        db_json, df_sdf=df_sdf, name_column="m_abbr", mol_colname="ROMol"
    )

    with open(out, "w") as fp:
        json.dump(db_json_augmented, fp)

    return db_json_augmented


if __name__ == "__main__":
    app()
