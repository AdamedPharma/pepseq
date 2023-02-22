import math
from io import BytesIO

import cairo

aa_color_dict = {
    'F': {
        'hexcolor': '#649DA0',
        'rgb_fractions': (0.390625, 0.61328125, 0.625)
        },
    'Y': {
        'hexcolor': '#7cbfb6',
        'rgb_fractions': (0.484375, 0.74609375, 0.7109375)
        },
    'W': {
        'hexcolor': '#5f9583',
        'rgb_fractions': (0.37109375, 0.58203125, 0.51171875)
        },
    'I': {
        'hexcolor': '#698B69',
        'rgb_fractions': (0.41015625, 0.54296875, 0.41015625)
        },
    'L': {
        'hexcolor': '#7daa7d',
        'rgb_fractions': (0.48828125, 0.6640625, 0.48828125)
        },
    'V': {
        'hexcolor': '#9BCD9B',
        'rgb_fractions': (0.60546875, 0.80078125, 0.60546875)
        },
    'A': {
        'hexcolor': '#a5d4a5',
        'rgb_fractions': (0.64453125, 0.828125, 0.64453125)
        },
    'X': {
        'hexcolor': '#de7f3b',
        'rgb_fractions': (0.8671875, 0.49609375, 0.23046875)
        },
    'M': {
        'hexcolor': '#80373b',
        'rgb_fractions': (0.5, 0.21484375, 0.23046875)
        },
    'C': {
        'hexcolor': '#95325b',
        'rgb_fractions': (0.58203125, 0.1953125, 0.35546875)
        },
    'P': {
        'hexcolor': '#e2b540',
        'rgb_fractions': (0.8828125, 0.70703125, 0.25)
        },
    'G': {
        'hexcolor': '#dcc750',
        'rgb_fractions': (0.859375, 0.77734375, 0.3125)
        },
    'T': {
        'hexcolor': '#808cbc',
        'rgb_fractions': (0.5, 0.546875, 0.734375)
        },
    'S': {
        'hexcolor': '#abb7d8',
        'rgb_fractions': (0.66796875, 0.71484375, 0.84375)
        },
    'Q': {
        'hexcolor': '#74badc',
        'rgb_fractions': (0.453125, 0.7265625, 0.859375)
        },
    'N': {
        'hexcolor': '#6ab7d2',
        'rgb_fractions': (0.4140625, 0.71484375, 0.8203125)
        },
    'H': {
        'hexcolor': '#664a73',
        'rgb_fractions': (0.3984375, 0.2890625, 0.44921875)
        },
    'R': {
        'hexcolor': '#104E8B',
        'rgb_fractions': (0.0625, 0.3046875, 0.54296875)
        },
    'K': {
        'hexcolor': '#1c407d',
        'rgb_fractions': (0.109375, 0.25, 0.48828125)
        },
    'E': {
        'hexcolor': '#c73333',
        'rgb_fractions': (0.77734375, 0.19921875, 0.19921875)
        },
    'D': {
        'hexcolor': '#d74242',
        'rgb_fractions': (0.83984375, 0.2578125, 0.2578125)
        },
    '-': {
        'hexcolor': '#bcbab3',
        'rgb_fractions': (0.734375, 0.7265625, 0.69921875)
        }
    }


def get_start_x(
    left_margin=100, is_corner=None, forward=None, step_x=None, is_start=None
):

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


def generate_kwargs_for_text_in_ellipse_balls(symbols, y, forward=True,
                                              is_corner=False, is_start=False, left_margin=100):
    step_x = 100
    startx = get_start_x(
                left_margin=left_margin,
                is_corner=is_corner,
                forward=forward,
                step_x=step_x,
                is_start=is_start)
    x = startx
    kwargs_list = []

    for symbol_position in range(len(symbols)):
        symbol = symbols[symbol_position]
        if len(symbol) == 1:
            font_size = 34
            text_x = x - 15
        else:
            font_size = 22
            text_x = x - 47
        kwargs = {
            'x': text_x,
            'y': y,
            'font_size': font_size,
            'text': symbol
            }
        kwargs_list.append(kwargs)
        if forward:
            x += step_x
        else:
            x -= step_x
    return kwargs_list


def generate_kwargs_for_ellipse_balls(symbols, y, forward=True, is_corner=False, is_start=False, left_margin=100):
    step_x = 100
    startx = get_start_x(
                left_margin=left_margin,
                is_corner=is_corner,
                forward=forward,
                step_x=step_x,
                is_start=is_start)
    x = startx
    kwargs_list = []

    for symbol_position in range(len(symbols)):
        symbol = symbols[symbol_position]
        symbol_dict = aa_color_dict.get(symbol)

        if symbol_dict is not None:
            rgb_fractions = symbol_dict['rgb_fractions']
            outline_rgb_fractions = (0.3, 0.3, 0.3)
            outline_width = 2
        else:
            rgb_fractions = (0.3, 0.3, 0.3)
            outline_rgb_fractions = (1.0, 0.0, 0.0)
            outline_width = 6
        kwargs = {
            'y': y,
            'x': x,
            'rgb_fractions': rgb_fractions,
            'outline_rgb_fractions': outline_rgb_fractions,
            'outline_width': outline_width
            }
        kwargs_list.append(kwargs)
        if forward:
            x += step_x
        else:
            x -= step_x
    return kwargs_list


def get_is_corner(num_iteration=None):
    """
    returns: corner or non_corner
    """
    """
    if num_iteration == 0:
        return False
    """
    odd = bool(num_iteration % 2)
    return odd


def get_direction(num_iteration=None):
    remainder = num_iteration % 4
    if remainder in [0, 3]:
        return "forward"
    elif remainder in [1, 2]:
        return "reverse"


def get_fragment_length(remaining_length=None, num_iteration=None):
    if num_iteration == 0:
        fragment_length = 9
    else:
        remainder = num_iteration % 2

        if remainder == 1:
            fragment_length = 1
        elif remainder == 0:
            fragment_length = 8
    return min(remaining_length, fragment_length)


def schema_layout_generator_from_symbols(symbols):

    remaining_symbols = symbols
    remaining_symbols.reverse()

    num_iteration = 0

    while True:
        remaining_length = len(remaining_symbols)

        if remaining_length <= 0:
            break

        fragment_direction = get_direction(num_iteration=num_iteration)
        fragment_length = get_fragment_length(
            remaining_length=remaining_length, num_iteration=num_iteration
        )
        seq_fragment = []

        for i in range(fragment_length):
            seq_fragment.append(remaining_symbols.pop())

        is_corner = get_is_corner(num_iteration=num_iteration)
        yield seq_fragment, fragment_length, fragment_direction, is_corner

        num_iteration += 1


def get_fragment_kwargs(symbols, y=70, is_start=True, fragment_direction='forward',
                        is_corner=False):
    if fragment_direction == 'forward':
        forward = True
    elif fragment_direction == 'reverse':
        forward = False

    kwargsy = generate_kwargs_for_ellipse_balls(symbols, y, forward, is_corner, is_start)
    kwargsy_text = generate_kwargs_for_text_in_ellipse_balls(symbols, y+10, forward,
                                                             is_corner, is_start)
    return kwargsy, kwargsy_text


def get_kwargs_from_symbols(symbols):
    fragments = list(schema_layout_generator_from_symbols(symbols))
    y = 70
    is_start = True

    all_kwargs_list = []
    all_kwargs_text_list = []

    for fragment in fragments:
        seq_fragment, fragment_length, fragment_direction, is_corner = fragment

        kwargs_list, text_kwargs_list = get_fragment_kwargs(seq_fragment, y, is_start=is_start,
                                                            fragment_direction=fragment_direction, is_corner=is_corner)
        all_kwargs_list += kwargs_list
        all_kwargs_text_list += text_kwargs_list
        is_start = False
        y += 55
    return all_kwargs_list, all_kwargs_text_list


def draw_ellipse_ball(cairo_context: cairo.Context, x: int, y: int, rgb_fractions: tuple,
                      outline_rgb_fractions=(0.3, 0.3, 0.3), outline_width=2):
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
    r = 50 + outline_width
    cairo_context.arc(0, 0, r, 0, 2 * math.pi)
    cairo_context.fill()
    cairo_context.set_source_rgb(*rgb_fractions)
    cairo_context.arc(0, 0, 50, 0, 2 * math.pi)
    cairo_context.fill()
    cairo_context.restore()
    return cairo_context


def draw_ellipse_balls(cairo_context: cairo.Context, keyword_args_sequence: list):
    """
    draws a chain of ellipses on cairo.Context provided
    a sequence of keyword argument dictionaries

        dict(x: int, y: int, rgb_fractions: tuple)

    The color of ellipse inside is given by rgb_fractions

    Outline of ellipse is grey (0.3, 0.3, 0.3)

    """
    for kwargs in keyword_args_sequence:
        x = kwargs['x']
        y = kwargs['y']
        rgb_fractions = kwargs['rgb_fractions']
        outline_rgb_fractions = kwargs.get('outline_rgb_fractions', (0.3, 0.3, 0.3))
        outline_width = kwargs.get('outline_width', 2)
        cairo_context = draw_ellipse_ball(cairo_context, x, y, rgb_fractions, outline_rgb_fractions, outline_width)
    return cairo_context


def draw_text_in_ellipse(cairo_context, x: int, y: int, text: str, font_size=34):
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


def draw_text_in_ellipse_balls(cairo_context: cairo.Context, keyword_args_sequence: list):
    """
    draws text inside a chain of ellipses on cairo.Context provided
    a sequence of keyword argument dictionaries

        dict(x: int, y: int, text: str, font_size=34)

    The color of text is white
    and font is Purisa

    """
    for kwargs in keyword_args_sequence:
        x = kwargs['x']
        y = kwargs['y']
        text = kwargs['text']
        font_size = kwargs['font_size']
        cairo_context = draw_text_in_ellipse(cairo_context, x, y, text, font_size)
    return cairo_context


def get_png_string_from_surface(surface):
    buffer = BytesIO()
    surface.write_to_png(buffer)  # _stream
    pngData = buffer.getvalue()
    buffer.close()
    return pngData


def draw_symbols(symbols, width=1024, height=1024):
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    cairo_context = cairo.Context(surface)
    kwargs_list, kwargs_text_list = get_kwargs_from_symbols(symbols)
    draw_ellipse_balls(cairo_context, kwargs_list)
    draw_text_in_ellipse_balls(cairo_context, kwargs_text_list)
    png_string = get_png_string_from_surface(surface)

    surface.flush()
    surface.finish()
    return png_string
