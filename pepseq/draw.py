from pepseq.Drawer import draw_symbols
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json



def draw_pepseq(
    pepseq: str,
    width: int = 1124,
    height: int = 640,
    omit_standard_termini: bool = False,
) -> str:
    """
    Draw a visual representation of a peptide sequence.
    Draw Peptide Schema Image based on sequence. Save it to (image_file_name)
    if image_file_name is not provided an UUID is generated.

    :param pepseq: The peptide sequence to be drawn.
    :type  pepseq: str
    :param width: The width of the drawing in pixels. Defaults to 1124.
    :type  width: int
    :param height: The height of the drawing in pixels. Defaults to 640.
    :type  height: int
    :param omit_standard_termini: Whether to omit standard termini (N-terminus and C-terminus) from the drawing. Defaults to False.
    :type  omit_standard_termini: bool

    :return: The PNG image of the drawn peptide sequence.
    :rtype: str
    """

    symbols = get_pep_json(pepseq)["symbols"]
    indices_to_exclude = []

    termini_present = set(["N", "C"])

    if omit_standard_termini:
        if symbols[0] == "H":
            indices_to_exclude.append(0)
            termini_present.remove("N")
        last_index = len(symbols) - 1
        if symbols[last_index] == "OH":
            indices_to_exclude.append(last_index)
            termini_present.remove("C")

    termini_present = list(termini_present)

    new_symbols = [
        symbols[i] for i in range(len(symbols)) if i not in indices_to_exclude
    ]

    png_string = draw_symbols(
        new_symbols, width=width, height=height, termini_present=termini_present
    )
    return png_string
