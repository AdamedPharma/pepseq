import math
import copy
from io import BytesIO
from collections import ChainMap

import cairo

aa_color_dict = {
    "F": {"hexcolor": "#649DA0", "rgb_fractions": (0.390625, 0.61328125, 0.625)},
    "Y": {"hexcolor": "#7cbfb6", "rgb_fractions": (0.484375, 0.74609375, 0.7109375)},
    "W": {"hexcolor": "#5f9583", "rgb_fractions": (0.37109375, 0.58203125, 0.51171875)},
    "I": {"hexcolor": "#698B69", "rgb_fractions": (0.41015625, 0.54296875, 0.41015625)},
    "L": {"hexcolor": "#7daa7d", "rgb_fractions": (0.48828125, 0.6640625, 0.48828125)},
    "V": {"hexcolor": "#9BCD9B", "rgb_fractions": (0.60546875, 0.80078125, 0.60546875)},
    "A": {"hexcolor": "#a5d4a5", "rgb_fractions": (0.64453125, 0.828125, 0.64453125)},
    "X": {"hexcolor": "#de7f3b", "rgb_fractions": (0.8671875, 0.49609375, 0.23046875)},
    "M": {"hexcolor": "#80373b", "rgb_fractions": (0.5, 0.21484375, 0.23046875)},
    "C": {"hexcolor": "#95325b", "rgb_fractions": (0.58203125, 0.1953125, 0.35546875)},
    "P": {"hexcolor": "#e2b540", "rgb_fractions": (0.8828125, 0.70703125, 0.25)},
    "G": {"hexcolor": "#dcc750", "rgb_fractions": (0.859375, 0.77734375, 0.3125)},
    "T": {"hexcolor": "#808cbc", "rgb_fractions": (0.5, 0.546875, 0.734375)},
    "S": {"hexcolor": "#abb7d8", "rgb_fractions": (0.66796875, 0.71484375, 0.84375)},
    "Q": {"hexcolor": "#74badc", "rgb_fractions": (0.453125, 0.7265625, 0.859375)},
    "N": {"hexcolor": "#6ab7d2", "rgb_fractions": (0.4140625, 0.71484375, 0.8203125)},
    "H": {"hexcolor": "#664a73", "rgb_fractions": (0.3984375, 0.2890625, 0.44921875)},
    "R": {"hexcolor": "#104E8B", "rgb_fractions": (0.0625, 0.3046875, 0.54296875)},
    "K": {"hexcolor": "#1c407d", "rgb_fractions": (0.109375, 0.25, 0.48828125)},
    "E": {"hexcolor": "#c73333", "rgb_fractions": (0.77734375, 0.19921875, 0.19921875)},
    "D": {"hexcolor": "#d74242", "rgb_fractions": (0.83984375, 0.2578125, 0.2578125)},
    "-": {"hexcolor": "#bcbab3", "rgb_fractions": (0.734375, 0.7265625, 0.69921875)},
}




def get_start_x(
    left_margin: int = 100, is_corner: bool = None, forward: bool = None, step_x: float = None,
    is_start: bool = None
) -> float:
    """
    Get starting X coordinate for 


     ACDEF
          G
     NMLKI
    P
     QRSTV
          W
         Y

    """

    right_margin = left_margin + (8.6 * step_x)
    last_right = left_margin + (8 * step_x)

    if forward:
        if is_corner:
            return left_margin + (0.4 * step_x)
        elif is_start:
            return left_margin
        else:
            return left_margin + (1.0 * step_x)
    else:
        if is_corner:
            return right_margin
        else:
            return last_right


def get_rev_x_and_font_size(symbol: str, variable_font_size: int = True) -> tuple:
    """
    """
    length = len(symbol)
    if not variable_font_size:
        if length == 1:
            font_size = 34
            rev_x = 15

        elif length == 3:
            font_size = 22
            rev_x = 30

        else:
            font_size = 22
            rev_x = 45
    else:
        f = (21/34)**(1/5)
        font_size = round(34 * (f**(length-1)))
        d_rev_x = {
            1 : (15+0),
            2: (15+7),
            3: (15+15),
            4: (15+18),
            5: (15+24),
            6: (15+32)
            }

        rev_x = d_rev_x.get(length, 47)

    return rev_x, font_size


def generate_kwargs_for_text_in_ellipse_balls(
    symbols: list,
    y: float,
    forward: bool = True,
    is_corner: bool = False,
    is_start: bool = False,
    left_margin: float = 100,
    step_x: float = 100
) -> list:
    """
    """
    
    startx = get_start_x(
        left_margin=left_margin,
        is_corner=is_corner,
        forward=forward,
        step_x=step_x,
        is_start=is_start,
    )
    x = startx
    kwargs_list = []

    

    for symbol_position in range(len(symbols)):
        symbol = symbols[symbol_position]

        rev_x, font_size = get_rev_x_and_font_size(symbol)
        text_x = x - rev_x
        if is_start and (symbol_position ==0):
            text_x += 15


        kwargs = {
            "x": text_x, "y": y, "font_size": font_size, "text": symbol}
        kwargs_list.append(kwargs)
        if forward:
            x += step_x
        else:
            x -= step_x

    return kwargs_list


def generate_kwargs_for_ellipse_balls(
    symbols: list, y: float, forward: bool = True, is_corner: bool = False,
    is_start: bool = False, left_margin: float = 100
) -> list:
    """
    """
    step_x = 100
    startx = get_start_x(
        left_margin=left_margin,
        is_corner=is_corner,
        forward=forward,
        step_x=step_x,
        is_start=is_start,
    )
    x = startx
    kwargs_list = []

    style = dict( radius=50 )

    modified_amino_acid_ellipse_style = dict(
            rgb_fractions = (0.3, 0.3, 0.3),
            outline_rgb_fractions = (1.0, 0.0, 0.0),
            outline_width = 6)
    
    loc = {}
    dict3 = dict(ChainMap(style, modified_amino_acid_ellipse_style))

    amino_acid_ellipse_style = dict(

            outline_rgb_fractions = (1.0, 0.0, 0.0),
            outline_width = 6)

    #set locations


    for symbol_position in range(len(symbols)):
        symbol = symbols[symbol_position]

        symbol_dict = aa_color_dict.get(symbol)

        if symbol_dict is not None:
            rgb_fractions = symbol_dict["rgb_fractions"]
            outline_rgb_fractions = (0.3, 0.3, 0.3)
            outline_width = 2
        else:
            rgb_fractions = (0.3, 0.3, 0.3)
            outline_rgb_fractions = (1.0, 0.0, 0.0)
            outline_width = 6
        kwargs = {
            "y": y,
            "x": x,
            "rgb_fractions": rgb_fractions,
            "outline_rgb_fractions": outline_rgb_fractions,
            "outline_width": outline_width,
            "radius": 50,
        }
        kwargs_list.append(kwargs)
        if forward:
            x += step_x
        else:
            x -= step_x
    return kwargs_list


def get_is_corner(num_iteration: int = None) -> bool:
    """
    ACDEFGHIK       num_iteration = 0 => False
             L      num_iteration = 1 => True
     VTSRQPNM       num_iteration = 2 => False
    W               num_iteration = 3 => True
     Y              num_iteration = 4 => False
    """

    odd = bool(num_iteration % 2)
    return odd


def get_direction(num_iteration: int = None) -> str:
    """
    ACDEFGHIK       num_iteration = 0 => forward
             L      num_iteration = 1 => reverse
     VTSRQPNM       num_iteration = 2 => reverse
    W               num_iteration = 3 => forward
     Y              num_iteration = 4 => forward
        
    """
    remainder = num_iteration % 4
    if remainder in [0, 3]:
        return "forward"
    elif remainder in [1, 2]:
        return "reverse"


def get_fragment_length(remaining_length: int=None, num_iteration: int=None,
                        lengths: dict={
                            'first': 9,
                            'corner': 1,
                            'standard': 8
                        }) -> int:
    """


    ACDEFGHIK       num_iteration = 0; remaining_length = 20 => length = 9
             L      num_iteration = 1; remaining_length = 11 => length = 1
     VTSRQPNM       num_iteration = 2; remaining_length = 10 => length = 8
    W               num_iteration = 3; remaining_length = 2 => length = 1
     Y              num_iteration = 4; remaining_length = 1 => length = 1
    """
    is_first = (num_iteration == 0)
    if is_first:
        
        fragment_length = lengths['first']
    else:
        odd_iteration = ((num_iteration % 2) == 1)
        is_corner = odd_iteration

        if is_corner:
            fragment_length = lengths['corner']
        else:
            fragment_length = lengths['standard']
    return min(remaining_length, fragment_length)


def schema_layout_generator_from_symbols(symbols: list):
    """
    input:

    symbols = [
        'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'
    ]
    
    is generator

    yields

    list(output) =   [
        (['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K'], 9, 'forward', False),
        (['L'], 1, 'reverse', True),
        (['M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V'], 8, 'reverse', False),
        (['W'], 1, 'forward', True),
        (['Y'], 1, 'forward', False)
        ]

    """
    remaining_symbols = symbols
    remaining_symbols.reverse()

    num_iteration = 0

    while True:
        remaining_length = len(remaining_symbols)

        if remaining_length <= 0:
            break

        fragment_direction = get_direction(num_iteration=num_iteration)
        fragment_length = get_fragment_length(
            remaining_length=remaining_length, num_iteration=num_iteration, lengths ={
                            'first': 9,
                            'corner': 1,
                            'standard': 8
                        }
        )
        seq_fragment = []

        for i in range(fragment_length):
            seq_fragment.append(remaining_symbols.pop())

        is_corner = get_is_corner(num_iteration=num_iteration)
        yield seq_fragment, fragment_length, fragment_direction, is_corner

        num_iteration += 1


def get_fragment_kwargs(
    symbols: list, y: float = 70, is_start: bool = True, fragment_direction: str = "forward",
    is_corner: bool = False
) -> tuple:
    """
            seq_fragment, fragment_length, fragment_direction, is_corner = fragment

    """
    if fragment_direction == "forward":
        forward = True
    elif fragment_direction == "reverse":
        forward = False

    kwargsy = generate_kwargs_for_ellipse_balls(
        symbols, y, forward, is_corner, is_start
    )
    kwargsy_text = generate_kwargs_for_text_in_ellipse_balls(
        symbols, y + 10, forward, is_corner, is_start
    )
    return kwargsy, kwargsy_text


def get_N_terminus_params(params: dict) -> dict:
    """
    gets N terminus params
    """
    out_params = copy.deepcopy(params)
    mod_params = {
        "radius": 35,
        "outline_width": 0,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3)
        }
    out_params.update(mod_params)
    out_params["x"] += 15

    return out_params


def get_C_terminus_params(params: dict, previous_params:dict) -> dict:
    """
    gets C terminus params
    """
    c_terminus_x = params["x"]
    out_params = copy.deepcopy(params)

    previous_x = previous_params["x"]

    if c_terminus_x > previous_x:
        out_params['x'] -= 15
    else:
        out_params['x'] += 15

    mod_params = {
        "radius": 35,
        "outline_width": 0,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3)
        }
    
    out_params.update(mod_params)

    return out_params

def get_kwargs_from_symbols(symbols: list, termini_present: list = ["N", "C"]) -> tuple:
    """
    symbols = [
        'CH3', 'Y', 'aMeAla', 'Q', 'G', 'T', 'F', 'T', 'S', 'D', 'Y', 'S', 'K', 'Y', 'L'] + \
    ['D', 'E', 'Cys(R1)', 'A', 'A', 'K', 'D', 'F', 'V', 'Cys(R2)', 'W', 'L', 'L', 'D', 'H', 'H'] + \
    ['P', 'S', 'S', 'G', 'Q', 'P', 'P', 'P', 'S', 'CN', 'C', '123456789', 'NH2']

    """
    fragments = list(schema_layout_generator_from_symbols(symbols))
    y = 70
    is_start = True

    all_kwargs_list = []
    all_kwargs_text_list = []

    for fragment in fragments:
        seq_fragment, fragment_length, fragment_direction, is_corner = fragment

        kwargs_list, text_kwargs_list = get_fragment_kwargs(
            seq_fragment,
            y,
            is_start=is_start,
            fragment_direction=fragment_direction,
            is_corner=is_corner,
        )
        all_kwargs_list += kwargs_list
        all_kwargs_text_list += text_kwargs_list
        is_start = False
        y += 55
    
    if "N" in termini_present:
        all_kwargs_list[0] = get_N_terminus_params(all_kwargs_list[0])
        offset = max(3 - len(all_kwargs_text_list[0]["text"]), 0)
        all_kwargs_text_list[0]["x"] += offset * 10
        all_kwargs_text_list[0]["font_size"] = 22

    if "C" in termini_present:
        c_terminus_params = get_C_terminus_params(params=all_kwargs_list[-1], previous_params=all_kwargs_list[-2])

        all_kwargs_list[-1] = c_terminus_params

        all_kwargs_text_list[-1]["x"] = all_kwargs_list[-1]["x"] - 30
        offset = max(3 - len(all_kwargs_text_list[-1]["text"]), 0)
        all_kwargs_text_list[-1]["x"] += offset * 10
        all_kwargs_text_list[-1]["font_size"] = 22

    middle_index = int(len(all_kwargs_list) / 2)
    middle_y = all_kwargs_list[middle_index]["y"]
    target_y = 320
    y_offset = target_y - middle_y

    for i in range(len(all_kwargs_text_list)):
        all_kwargs_text_list[i]["y"] += y_offset
        all_kwargs_list[i]["y"] += y_offset

    return all_kwargs_list, all_kwargs_text_list


def draw_ellipse_ball(
    cairo_context: cairo.Context,
    x: int,
    y: int,
    rgb_fractions: tuple,
    outline_rgb_fractions: tuple = (0.3, 0.3, 0.3),
    outline_width: int = 2,
    radius: int = 50,
) -> cairo.Context:
    """
    draws an ellipse on cairo.Context provided
    at x and y coordinates

    The color of ellipse inside is given by rgb_fractions

    Outline of ellipse is grey (0.3, 0.3, 0.3)

    """

    cairo_context.save()
    cairo_context.translate(x, y)
    cairo_context.scale(1, 0.7)
    cairo_context.set_source_rgb(*outline_rgb_fractions)
    r = radius + outline_width
    cairo_context.arc(0, 0, r, 0, 2 * math.pi)
    cairo_context.fill()
    cairo_context.set_source_rgb(*rgb_fractions)
    cairo_context.arc(0, 0, radius, 0, 2 * math.pi)
    cairo_context.fill()
    cairo_context.restore()
    return cairo_context


def draw_ellipse_balls(cairo_context: cairo.Context, keyword_args_sequence: list
                       ) -> cairo.Context:
    """
    draws a chain of ellipses on cairo.Context provided
    a sequence of keyword argument dictionaries

        dict(x: int, y: int, rgb_fractions: tuple)

    The color of ellipse inside is given by rgb_fractions

    Outline of ellipse is grey (0.3, 0.3, 0.3)

    """
    for kwargs in keyword_args_sequence:
        x = kwargs["x"]
        y = kwargs["y"]
        radius = kwargs["radius"]
        rgb_fractions = kwargs["rgb_fractions"]
        outline_rgb_fractions = kwargs.get("outline_rgb_fractions", (0.3, 0.3, 0.3))
        outline_width = kwargs.get("outline_width", 2)
        cairo_context = draw_ellipse_ball(
            cairo_context,
            x,
            y,
            rgb_fractions,
            outline_rgb_fractions,
            outline_width,
            radius=radius,
        )
    return cairo_context


def draw_text_in_ellipse(cairo_context: cairo.Context, x: int, y: int, text: str,
                          font_size=34) -> cairo.Context:
    """
    writes text in ellipse on cairo.Context
    at x and y coordinates

    The color of text is white
    and font is Purisa

    """
    text_rgb_fractions = (1.0, 1.0, 1.0)
    cairo_context.set_source_rgb(*text_rgb_fractions)
    font_face_args = ("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    cairo_context.select_font_face(*font_face_args)
    cairo_context.set_font_size(font_size)
    cairo_context.move_to(x, y)
    cairo_context.show_text(text)
    return cairo_context


def draw_text_in_ellipse_balls(
    cairo_context: cairo.Context, keyword_args_sequence: list
) -> cairo.Context:
    """
    draws text inside a chain of ellipses on cairo.Context provided
    a sequence of keyword argument dictionaries

        dict(x: int, y: int, text: str, font_size=34)

    The color of text is white
    and font is Purisa

    """
    for kwargs in keyword_args_sequence:
        x = kwargs["x"]
        y = kwargs["y"]
        text = kwargs["text"]
        font_size = kwargs["font_size"]
        cairo_context = draw_text_in_ellipse(cairo_context, x, y, text, font_size)
    return cairo_context


def get_png_string_from_surface(surface: cairo.ImageSurface) -> bytes:
    """
    turns cairo.ImageSurface into png byte buffer
    and then into PNG image format string
    """
    buffer = BytesIO()
    surface.write_to_png(buffer)  # _stream
    pngData = buffer.getvalue()
    buffer.close()
    return pngData


def draw_symbols(symbols: list, width: int = 1024, height: int = 1024,
                  termini_present: list = ["N", "C"], out: str = None) -> str:
    """
    From symbols list (a sequence already split into residue symbols)
    Input:
        symbols: residue symbols

        symbols = [
             'CH3', 'Y', 'aMeAla', 'Q', 'G', 'T', 'F', 'T', 'S', 'D', 'Y', 'S', 'K', 'Y', 'L'] + \
            ['D', 'E', 'Cys(R1)', 'A', 'A', 'K', 'D', 'F', 'V', 'Cys(R2)', 'W', 'L', 'L', 'D', 'H', 'H'] + \
            ['P', 'S', 'S', 'G', 'Q', 'P', 'P', 'P', 'S', 'CN', 'C', 'R', 'NH2']
        
        width: image width (px)
        height: image height (px)
        termini_present:
            whether N terminus different than 'H' is present
            and/or 
            whether C terminus different than 'OH' is present

    """
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    cairo_context = cairo.Context(surface)
    kwargs_list, kwargs_text_list = get_kwargs_from_symbols(
        symbols, termini_present=termini_present
    )
    draw_ellipse_balls(cairo_context, kwargs_list)
    draw_text_in_ellipse_balls(cairo_context, kwargs_text_list)

    if out is not None:
        surface.write_to_png(out)
        surface.flush()
        surface.finish()
        return out
    else:

        png_string = get_png_string_from_surface(surface)

        surface.flush()
        surface.finish()
        return png_string
