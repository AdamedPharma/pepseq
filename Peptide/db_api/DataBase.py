import json
from collections import namedtuple

AminoAcid = namedtuple(
    "AminoAcid", "smiles smiles_radical three_letter_code smarts smarts_comments"
)


def read_AminoAcid_from_dict(in_dict):

    """
    creates and returns AminoAcid namedtuple object
    from dictionary object provided
    """

    return AminoAcid(
        in_dict["smiles"],
        in_dict["smiles_radical"],
        in_dict["three_letter_code"],
        in_dict["smarts"],
        in_dict["smarts_comments"],
    )


def read_aa_smiles_dict(in_dict):
    """
    creates and returns a dictionary of AminoAcid namedtuple objects
    assigned to amino acid name
        from dictionary object provided
    """

    out_dict = {}

    for amino_acid_name in in_dict:
        amino_acid_dict = in_dict[amino_acid_name]
        out_dict[amino_acid_name] = read_AminoAcid_from_dict(amino_acid_dict)

    return out_dict


def lower_dict(in_dict=None):
    """

    transforms dictionary:
        {
            "Ala": "A",
            "Cys": "C",
            }
    to:
        {
            "ala": "a",
            "cys": "c",
            }

    """

    out_dict = {}
    for k in in_dict:
        val = in_dict[k]
        k_lower = k.lower()
        out_dict[k_lower] = val.lower()
    return out_dict


def join_dicts(dicts):
    """

    dicts:

    [   {
            "Ala": "A",
            "Cys": "C",
            },
        {
            "ala": "a",
            "cys": "c",
            },
        {
            "dA": "a",
            "dC": "c",
            "dD": "d",
            "dE": "e",
            "dAla": "a",
            "dCys": "c",
            "dAsp": "d"
            },
        {
            "aMeAla": "X",
            "2-Aminoisobutyric acid": "X",
            "2-Methylalanine": "X"
            }
        ]

    are joined into:

    {
        "Ala": "A",
        "Cys": "C",
        "ala": "a",
        "cys": "c",
        "dA": "a",
        "dC": "c",
        "dD": "d",
        "dE": "e",
        "dAla": "a",
        "dCys": "c",
        "dAsp": "d"
        "aMeAla": "X",
        "2-Aminoisobutyric acid": "X",
        "2-Methylalanine": "X"
        }
    """

    joined_dict = {}

    for dict2 in dicts:
        joined_dict.update(dict2)

    return joined_dict


class FileSystemDbRepo(object):
    """
    this represents a Facade for Reading From a JSON File
    condaining information on Amino Acid Species:
        their SMILES and SMARTS codes, etc.
    """

    def __init__(self):
        return

    @classmethod
    def read_from_json(cls, path=None, **kwargs):
        with open(path) as fp:
            db_json = json.load(fp)
        obj = cls()
        obj.read_json(db_json=db_json)
        obj.process_json()
        return obj

    def read_json(self, db_json=None):
        """
        reads in db_json
        """
        self.n_terms_smi_codes = db_json["n_terms_smi_codes"]
        self.n_terms_order = db_json["n_terms_order"]
        self.c_terms_smi_codes = db_json["c_terms_smi_codes"]
        self.c_terms_order = db_json["c_terms_order"]
        self.substitute_smiles = db_json["substitute_smiles"]
        self.substitute_dict = db_json["substitute_dict"]
        self.aa_order = db_json["aa_order"]
        self.mod_smi_codes = db_json.get("mod_smi_codes")

        self.l_proteogenic_3letter = db_json["l_proteogenic_3letter"]

        self.d_proteogenic_2letter = db_json["d_proteogenic_2letter"]
        self.d_proteogenic_4letter = db_json["d_proteogenic_4letter"]

        self.modified_aa_codes = db_json["modified_aa_codes"]

        self.file_aa_smiles_dict = db_json["aa_smiles_dict"]

    def process_json(self):

        self.l_proteogenic = set(self.file_aa_smiles_dict.keys())
        self.d_proteogenic = set([i.lower() for i in self.l_proteogenic])

        self.aa_smiles_dict = read_aa_smiles_dict(self.file_aa_smiles_dict)

        self.d_proteogenic_3letter = lower_dict(in_dict=self.l_proteogenic_3letter)

        self.three_letter = join_dicts(
            [self.l_proteogenic_3letter, self.d_proteogenic_3letter]
        )
        aa_symbols = {}
        for key in self.aa_smiles_dict:
            pass

        self.symbols = join_dicts(
            [
                self.d_proteogenic_2letter,
                self.d_proteogenic_4letter,
                self.three_letter,
                self.modified_aa_codes,
            ]
        )

        return
