import math
import copy
from io import BytesIO

import cairo
from typing import Union


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
    left_margin: int = 100,
    is_corner: Union[bool, None] = None,
    forward: Union[bool, None] = None,
    step_x: Union[float, None] = None,
    is_start: Union[bool, None] = None,
) -> float:
    """
    Get starting X coordinate for next line of drawn sequence

    :param left_margin: int = left margin of the image
    :type left_margin: int

    :param is_corner: bool = whether the sequence is a corner
    :type is_corner: bool

    :param forward: bool = whether the sequence is drawn forward
    :type forward: bool

    :param step_x: float = step in X direction
    :type step_x: float

    :param is_start: bool = whether the sequence is the first in the line
    :type is_start: bool

    :return: X coordinate
    :rtype: float
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


def get_rev_x_and_font_size(symbol: str, variable_font_size: Union[int, bool] = True) -> tuple:
    """
    Get reverse X coordinate and font size for the symbol

    :param symbol : symbol
    :type symbol : str

    :param variable_font_size : whether font size is variable
    :type variable_font_size : bool

    :return: reverse X coordinate and font size
    :rtype: tuple
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
        f = (21 / 34) ** (1 / 5)
        font_size = round(34 * (f ** (length - 1)))
        d_rev_x = {
            1: (15 + 0),
            2: (15 + 7),
            3: (15 + 15),
            4: (15 + 18),
            5: (15 + 24),
            6: (15 + 32),
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
    step_x: float = 100,
) -> list:
    """
    Generate kwargs for text in ellipse balls

    :param symbols: list of symbols
    :type symbols: list

    :param y: y coordinate
    :type y: float

    :param forward: whether the sequence is drawn forward
    :type forward: bool

    :param is_corner: whether the sequence is a corner
    :type is_corner: bool

    :param is_start: whether the sequence is the first in the line
    :type is_start: bool

    :param left_margin: left margin of the image
    :type left_margin: float

    :param step_x: step in X direction
    :type step_x: float

    :return: list: list of kwargs
    :rtype: list
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
        if is_start and (symbol_position == 0):
            text_x += 15

        kwargs = {"x": text_x, "y": y, "font_size": font_size, "text": symbol}
        kwargs_list.append(kwargs)
        if forward:
            x += step_x
        else:
            x -= step_x

    return kwargs_list


def generate_kwargs_for_ellipse_balls(
    symbols: list,
    y: float,
    forward: bool = True,
    is_corner: bool = False,
    is_start: bool = False,
    left_margin: float = 100,
) -> list:
    """
    Generate kwargs for ellipse balls

    :param symbols: list of symbols
    :type symbols: list

    :param y: y coordinate
    :type y: float

    :param forward: whether the sequence is drawn forward
    :type forward: bool

    :param is_corner: whether the sequence is a corner
    :type is_corner: bool

    :param is_start: whether the sequence is the first in the line
    :type is_start: bool

    :param left_margin: left margin of the image
    :type left_margin: float

    :return: list: list of kwargs
    :rtype: list
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

    # set locations

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


def get_is_corner(num_iteration: Union[int, None] = None) -> bool:
    """
    Every odd iteration is a corner and even is not corner.

    :param: num_iteration
    :type: int

    :return: whether the iteration is a corner
    :rtype: bool
    """

    odd = bool(num_iteration % 2)
    return odd


def get_direction(num_iteration: Union[int, None] = None) -> str:
    """
    Every 4th iteration is forward, 1st is reverse, 2nd is reverse, 3rd is forward

    :param num_iteration
    :type num_iteration

    :return: direction of the iteration
    :rtype: str
    """
    remainder = num_iteration % 4
    if remainder in [0, 3]:
        return "forward"
    elif remainder in [1, 2]:
        return "reverse"


def get_fragment_length(
    remaining_length: Union[int, None] = None,
    num_iteration: Union[int, None] = None,
    lengths: dict = {"first": 9, "corner": 1, "standard": 8},
) -> int:
    """
    Get length of sequence fragment. First fragment is 9, corner is 1, standard is 8.
    If length of remaining sequence is less than fragment length,
      fragment length is the length of the remaining sequence.

    :param remaining_length: length of the remaining sequence
    :type remaining_length: int

    :param num_iteration: iteration number
    :type num_iteration: int

    :param lengths: dictionary of lengths for first, corner and standard fragments
    :type lengths: dict

    :return: length of the fragment
    :rtype: int
    """
    is_first = num_iteration == 0
    if is_first:
        fragment_length = lengths["first"]
    else:
        odd_iteration = (num_iteration % 2) == 1
        is_corner = odd_iteration

        if is_corner:
            fragment_length = lengths["corner"]
        else:
            fragment_length = lengths["standard"]
    return min(remaining_length, fragment_length)


def schema_layout_generator_from_symbols(symbols: list):
    """
    Generate schema layout based on list of symbols. Schema layout is a list of tuples.
    Each tuple contains a list of symbols, length of the fragment,
      direction of the fragment and whether the fragment is a corner.

    :param symbols: list of symbols
    :type symbols: list

    :return: generator: generator of schema layout
    :rtype: generator
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
            remaining_length=remaining_length,
            num_iteration=num_iteration,
            lengths={"first": 9, "corner": 1, "standard": 8},
        )
        seq_fragment = []

        for i in range(fragment_length):
            seq_fragment.append(remaining_symbols.pop())

        is_corner = get_is_corner(num_iteration=num_iteration)
        yield seq_fragment, fragment_length, fragment_direction, is_corner

        num_iteration += 1


def get_fragment_kwargs(
    symbols: list,
    y: float = 70,
    is_start: bool = True,
    fragment_direction: str = "forward",
    is_corner: bool = False,
) -> tuple:
    """
    :param symbols: list of symbols
    :type symbols: list

    :param y: y coordinate
    :type y: float

    :param is_start: whether the sequence is the first in the line
    :type is_start: bool

    :param fragment_direction: direction of the fragment
    :type fragment_direction: str

    :param is_corner: whether the fragment is a corner
    :type is_corner: bool

    :return: tuple = tuple of kwargs for ellipse balls and text in ellipse balls
    seq_fragment, fragment_length, fragment_direction, is_corner = fragment
    :rtype: tuple
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

    :param params: dictionary of parameters
    :type params: dict

    :return: dictionary of parameters
    :rtype: dict
    """
    out_params = copy.deepcopy(params)
    mod_params = {
        "radius": 35,
        "outline_width": 0,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
    }
    out_params.update(mod_params)
    out_params["x"] += 15

    return out_params


def get_C_terminus_params(params: dict, previous_params: dict) -> dict:
    """
    gets C terminus params

    :param params: dictionary of parameters
    :type params: dict

    :param previous_params: dictionary of previous parameters
    :type previous_params: dict

    :return: dictionary of parameters
    :rtype: dict
    """

    c_terminus_x = params["x"]
    out_params = copy.deepcopy(params)

    previous_x = previous_params["x"]

    if c_terminus_x > previous_x:
        out_params["x"] -= 15
    else:
        out_params["x"] += 15

    mod_params = {
        "radius": 35,
        "outline_width": 0,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
    }

    out_params.update(mod_params)

    return out_params


def get_kwargs_from_symbols(symbols: list, termini_present: list = ["N", "C"]) -> tuple:
    """
    Get kwargs from symbols

    :param symbols: list of symbols
    :type symbols: list

    :param termini_present: whether N terminus different than 'H' and whether C terminus different than 'OH' is present
    :type termini_present: list

    :return: all_kwargs_list, all_kwargs_text_list: tuple of kwargs for ellipse balls and text in ellipse balls
    :rtype: tuple

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
        c_terminus_params = get_C_terminus_params(
            params=all_kwargs_list[-1], previous_params=all_kwargs_list[-2]
        )

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

    :param cairo_context: cairo context
    :type cairo_context: cairo.Context

    :param x: x coordinate
    :type x: int

    :param y: y coordinate
    :type y: int

    :param rgb_fractions: RGB fractions
    :type rgb_fractions: tuple

    :param outline_rgb_fractions: RGB fractions for outline
    :type outline_rgb_fractions: tuple

    :param outline_width: width of the outline
    :type outline_width: int

    :param radius: radius of the ellipse
    :type radius: int

    :return: cairo.Context: cairo context
    :rtype: cairo.Context
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


def draw_ellipse_balls(
    cairo_context: cairo.Context, keyword_args_sequence: list
) -> cairo.Context:
    """
    draws a chain of ellipses on cairo.Context provided
    a sequence of keyword argument dictionaries

        dict(x: int, y: int, rgb_fractions: tuple)

    The color of ellipse inside is given by rgb_fractions

    Outline of ellipse is grey (0.3, 0.3, 0.3)

    :param cairo_context: cairo context
    :type cairo_context: cairo.Context

    :param keyword_args_sequence: list of keyword argument dictionaries
    :type keyword_args_sequence: list

    :return: cairo context
    :rtype: cairo.Context
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


def draw_text_in_ellipse(
    cairo_context: cairo.Context, x: int, y: int, text: str, font_size=34
) -> cairo.Context:
    """
    writes text in ellipse on cairo.Context
    at x and y coordinates

    The color of text is white
    and font is Purisa

    :param cairo_context: cairo context
    :type cairo_context: cairo.Context

    :param x: x coordinate
    :type x: int

    :param y: y coordinate
    :type y: int

    :param text: text to be written
    :type text: str

    :param font_size: font size
    :type font_size: int

    :return: cairo context
    :rtype: cairo.Context

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

    :param cairo_context: cairo context
    :type cairo_context: cairo.Context

    :param keyword_args_sequence: list of keyword argument dictionaries
    :type keyword_args_sequence: list

    :return: cairo context
    :rtype: cairo.Context
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

    :param surface: cairo image surface
    :type surface: cairo.ImageSurface

    :return: PNG image format string
    :rtype: bytes
    """
    buffer = BytesIO()
    surface.write_to_png(buffer)  # _stream
    pngData = buffer.getvalue()
    buffer.close()
    return pngData


def draw_symbols(
    symbols: list,
    width: int = 1024,
    height: int = 1024,
    termini_present: list = ["N", "C"],
    out: Union[str, None] = None,
) -> str:
    """
    Example of symbols list is [
    'CH3', 'Y', 'aMeAla', 'Q', 'G', 'T', 'F', 'T', 'S', 'D', 'Y', 'S', 'K', 'Y', 'L'] + \
    ['D', 'E', 'Cys(R1)', 'A', 'A', 'K', 'D', 'F', 'V', 'Cys(R2)', 'W', 'L', 'L', 'D', 'H', 'H'] + \
    ['P', 'S', 'S', 'G', 'Q', 'P', 'P', 'P', 'S', 'CN', 'C', 'R', 'NH2']

    From symbols list (a sequence already split into residue symbols)

    :param symbols: residue symbols
    :type symbols: list

    :param width: image width (px)
    :type width: int

    :param height: image height (px)
    :type height: int

    :param termini_present: whether N terminus different than 'H'
      is present and/or whether C terminus different than 'OH' is present
    :type termini_present: list

    :param out: output file name
    :type out: str

    :return: PNG image format string
    :rtype: str
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
