from Peptide.utils.Parser import find_termini, parse_canonical2


def get_attachment_point_json(res_id, decomposition):
    ResName, attachment_point_id = decomposition
    d_atom_name = {"Cys": "SG", "Lys": "NZ"}

    if ResName in d_atom_name:
        AtomName = d_atom_name[ResName]
    else:
        AtomName = ""

    att_point_json = {
        "attachment_point_id": attachment_point_id,
        "ResID": str(res_id + 1),
        "AtomName": AtomName,
        "ResidueName": ResName,
    }
    return att_point_json


def decompose_symbol(s):
    if ("(" in s) and (")" in s):
        inside_bracket = s.split("(")[1].split(")")[0]
        before_bracket = s.split("(")[0]
        if "R" in inside_bracket:
            R_id = inside_bracket[1:]
            return before_bracket, R_id
    return s


def get_attachment_points_on_sequence_json(symbols):
    att_points = {}

    for res_id in range(len(symbols)):
        symbol = symbols[res_id]
        decomposition = decompose_symbol(symbol)
        if type(decomposition) == tuple:
            res_name, attachment_point_id = decomposition
            attachment_point_json = get_attachment_point_json(res_id, decomposition)
            att_points[int(attachment_point_id)] = attachment_point_json
    return att_points


def get_base_seq(symbols):
    three_to_one = {"Cys": "C", "Lys": "K"}

    base_seq = ""

    for res_id in range(len(symbols)):
        symbol = symbols[res_id]
        decomposition = decompose_symbol(symbol)
        if type(decomposition) == tuple:
            res_name, attachment_point_id = decomposition
            if res_name in three_to_one:
                res_name = three_to_one[res_name]
            symbol = res_name

        if len(symbol) > 1:
            symbol = "{%s}" % (symbol)
        base_seq = base_seq + symbol
    return base_seq


def get_ext_mod(pepseq, smiles="[1*]C[2*] |$R1;placeholder;R2$|"):
    symbols = parse_canonical2(pepseq)
    attachment_points_on_sequence = get_attachment_points_on_sequence_json(symbols)
    if attachment_points_on_sequence.keys():
        max_attachment_point_id = max(
            [int(i) for i in attachment_points_on_sequence.keys()]
        )

        ext_mod = {
            "smiles": smiles,
            "max_attachment_point_id": max_attachment_point_id,
            "attachment_points_on_sequence": attachment_points_on_sequence,
        }
        return ext_mod


def get_pep_json(pepseq_format, db_json, mod_smiles="[1*]C[2*] |$R1;placeholder;R2$|"):
    # N_terminus, pepseq, C_terminus = pepseq_format.split("~")
    N_terminus, C_terminus, pepseq = find_termini(pepseq_format, db_json)
    symbols = parse_canonical2(pepseq)
    base_seq = get_base_seq(symbols)
    ext_mod = get_ext_mod(pepseq, mod_smiles)

    pep_json = {
        "sequence": base_seq,
        "internal_modifications": [],
        "C_terminus": C_terminus,
        "N_terminus": N_terminus,
        "pepseq_format": pepseq_format,
    }
    if ext_mod is not None:
        pep_json["external_modifications"] = [ext_mod]
    else:
        pep_json["external_modifications"] = []

    return pep_json
