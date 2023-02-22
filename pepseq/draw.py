import os


from pepseq.Drawer import draw_symbols
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json


def draw_pepseq(
    pepseq: str,
    image_file_name: str = "out.png",
    width: int = 1124,
    height: int = 640,
):
    """
    Draw Peptide Schema Image based on sequence. Save it to (image_file_name)
    if image_file_name is not provided an UUID is generated.
    """

    symbols = get_pep_json(pepseq)["symbols"]
    symbols[0] = '%s-' % (symbols[0])
    symbols[-1] = '-%s' % (symbols[-1])
    draw_symbols(symbols, width=width, height=height, out=image_file_name)
    return


def draw_all(pepseq, root_dir='.', width: int = 1124,
             height: int = 640, image_file_name: str = "out.png"):
    os.system('mkdir -p %s' % root_dir)
    pepseq_filename = "%s/%s" % (root_dir, image_file_name)

    draw_pepseq(pepseq, pepseq_filename, width, height)
    return
