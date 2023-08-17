import pytest
from pepseq.Drawer import draw_symbols


def test_draw_sequence():
    symbols = ['Cys(R1)', 'S', 'Cys(R2)', 'S', 'D', 'E', 'K', 'S', 'D', 'E', 'K', 'M', 'N', 'P', 'Q'] + \
    ['Cys(R1)', 'S', 'Cys(R2)', 'S', 'D', 'E', 'K', 'S', 'D', 'E', 'K', 'M', 'N', 'P', 'Q', 'R']
    png_string = draw_symbols(symbols, width = 1024, height = 1024)#, out = '../image_file_name_kwargsy4.png')

    symbols = ['CH3', 'Y', 'aMeAla', 'Q', 'G', 'T', 'F', 'T', 'S', 'D', 'Y', 'S', 'K', 'Y', 'L'] + \
    ['D', 'E', 'Cys(R1)', 'A', 'A', 'K', 'D', 'F', 'V', 'Cys(R2)', 'W', 'L', 'L', 'D', 'H', 'H'] + \
    ['P', 'S', 'S', 'G', 'Q', 'P', 'P', 'P', 'S', 'NH2']
    png_string = draw_symbols(symbols, width = 1024, height = 1024)#, out = '../image_file_name_can_err.png')


    symbols = ['CH3', 'Y', 'aMeAla', 'Q', 'G', 'T', 'F', 'T', 'S', 'D', 'Y', 'S', 'K', 'Y', 'L'] + \
    ['D', 'E', 'Cys(R1)', 'A', 'A', 'K', 'D', 'F', 'V', 'Cys(R2)', 'W', 'L', 'L', 'D', 'H', 'H'] + \
    ['P', 'S', 'S', 'G', 'Q', 'P', 'P', 'P', 'S', 'CN', 'C', '123456789', 'NH2']
    png_string = draw_symbols(symbols, width = 1024, height = 1024)#, out = '../image_file_name_can_err_2.png')    

    symbols = ['CH3', 'Y', 'aMeAla', 'Q', 'G', 'T', 'F', 'T', 'S', 'D', 'Y', 'S', 'K', 'Y', 'L'] + \
    ['D', 'E', 'Cys(R1)', 'A', 'A', 'K', 'D', 'F', 'V', 'Cys(R2)', 'W', 'L', 'L', 'D', 'H', 'H'] + \
    ['P', 'S', 'S', 'G', 'Q', 'P', 'P', 'P', 'S', 'CN', 'C', '123456789', 'NH2']

    png_out = draw_symbols(symbols, width = 1024, height = 1024, out = '../image_file_name_can_err_2.png')    

    return