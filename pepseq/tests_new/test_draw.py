import copy

from pepseq.Drawer import (
    draw_symbols,
    get_kwargs_from_symbols,
    get_fragment_kwargs,
    generate_kwargs_for_text_in_ellipse_balls,
    generate_kwargs_for_ellipse_balls,
    get_rev_x_and_font_size,
    get_start_x,
    get_N_terminus_params,
    get_C_terminus_params,
)


problematic_symbols = (
    ["CH3", "Y", "aMeAla", "Q", "G", "T", "F", "T", "S", "D", "Y", "S", "K", "Y", "L"]
    + [
        "D",
        "E",
        "Cys(R1)",
        "A",
        "A",
        "K",
        "D",
        "F",
        "V",
        "Cys(R2)",
        "W",
        "L",
        "L",
        "D",
        "H",
        "H",
    ]
    + ["P", "S", "S", "G", "Q", "P", "P", "P", "S", "CN", "C", "123456789", "NH2"]
)


kwargs_list_as_tuple_mock = (
    (100, 115, (0.3, 0.3, 0.3), (0.3, 0.3, 0.3), 0, 35),
    (100, 200, (0.484375, 0.74609375, 0.7109375), (0.3, 0.3, 0.3), 2, 50),
    (100, 300, (0.3, 0.3, 0.3), (1.0, 0.0, 0.0), 6, 50),
    (100, 400, (0.453125, 0.7265625, 0.859375), (0.3, 0.3, 0.3), 2, 50),
    (100, 500, (0.859375, 0.77734375, 0.3125), (0.3, 0.3, 0.3), 2, 50),
    (100, 600, (0.5, 0.546875, 0.734375), (0.3, 0.3, 0.3), 2, 50),
    (100, 700, (0.390625, 0.61328125, 0.625), (0.3, 0.3, 0.3), 2, 50),
    (100, 800, (0.5, 0.546875, 0.734375), (0.3, 0.3, 0.3), 2, 50),
    (100, 900, (0.66796875, 0.71484375, 0.84375), (0.3, 0.3, 0.3), 2, 50),
    (155, 960.0, (0.83984375, 0.2578125, 0.2578125), (0.3, 0.3, 0.3), 2, 50),
    (210, 900, (0.484375, 0.74609375, 0.7109375), (0.3, 0.3, 0.3), 2, 50),
    (210, 800, (0.66796875, 0.71484375, 0.84375), (0.3, 0.3, 0.3), 2, 50),
    (210, 700, (0.109375, 0.25, 0.48828125), (0.3, 0.3, 0.3), 2, 50),
    (210, 600, (0.484375, 0.74609375, 0.7109375), (0.3, 0.3, 0.3), 2, 50),
    (210, 500, (0.48828125, 0.6640625, 0.48828125), (0.3, 0.3, 0.3), 2, 50),
    (210, 400, (0.83984375, 0.2578125, 0.2578125), (0.3, 0.3, 0.3), 2, 50),
    (210, 300, (0.77734375, 0.19921875, 0.19921875), (0.3, 0.3, 0.3), 2, 50),
    (210, 200, (0.3, 0.3, 0.3), (1.0, 0.0, 0.0), 6, 50),
    (265, 140.0, (0.64453125, 0.828125, 0.64453125), (0.3, 0.3, 0.3), 2, 50),
    (320, 200.0, (0.64453125, 0.828125, 0.64453125), (0.3, 0.3, 0.3), 2, 50),
    (320, 300.0, (0.109375, 0.25, 0.48828125), (0.3, 0.3, 0.3), 2, 50),
    (320, 400.0, (0.83984375, 0.2578125, 0.2578125), (0.3, 0.3, 0.3), 2, 50),
    (320, 500.0, (0.390625, 0.61328125, 0.625), (0.3, 0.3, 0.3), 2, 50),
    (320, 600.0, (0.60546875, 0.80078125, 0.60546875), (0.3, 0.3, 0.3), 2, 50),
    (320, 700.0, (0.3, 0.3, 0.3), (1.0, 0.0, 0.0), 6, 50),
    (320, 800.0, (0.37109375, 0.58203125, 0.51171875), (0.3, 0.3, 0.3), 2, 50),
    (320, 900.0, (0.48828125, 0.6640625, 0.48828125), (0.3, 0.3, 0.3), 2, 50),
    (375, 960.0, (0.48828125, 0.6640625, 0.48828125), (0.3, 0.3, 0.3), 2, 50),
    (430, 900, (0.83984375, 0.2578125, 0.2578125), (0.3, 0.3, 0.3), 2, 50),
    (430, 800, (0.3984375, 0.2890625, 0.44921875), (0.3, 0.3, 0.3), 2, 50),
    (430, 700, (0.3984375, 0.2890625, 0.44921875), (0.3, 0.3, 0.3), 2, 50),
    (430, 600, (0.8828125, 0.70703125, 0.25), (0.3, 0.3, 0.3), 2, 50),
    (430, 500, (0.66796875, 0.71484375, 0.84375), (0.3, 0.3, 0.3), 2, 50),
    (430, 400, (0.66796875, 0.71484375, 0.84375), (0.3, 0.3, 0.3), 2, 50),
    (430, 300, (0.859375, 0.77734375, 0.3125), (0.3, 0.3, 0.3), 2, 50),
    (430, 200, (0.453125, 0.7265625, 0.859375), (0.3, 0.3, 0.3), 2, 50),
    (485, 140.0, (0.8828125, 0.70703125, 0.25), (0.3, 0.3, 0.3), 2, 50),
    (540, 200.0, (0.8828125, 0.70703125, 0.25), (0.3, 0.3, 0.3), 2, 50),
    (540, 300.0, (0.8828125, 0.70703125, 0.25), (0.3, 0.3, 0.3), 2, 50),
    (540, 400.0, (0.66796875, 0.71484375, 0.84375), (0.3, 0.3, 0.3), 2, 50),
    (540, 500.0, (0.3, 0.3, 0.3), (1.0, 0.0, 0.0), 6, 50),
    (540, 600.0, (0.58203125, 0.1953125, 0.35546875), (0.3, 0.3, 0.3), 2, 50),
    (540, 700.0, (0.3, 0.3, 0.3), (1.0, 0.0, 0.0), 6, 50),
    (540, 785.0, (0.3, 0.3, 0.3), (0.3, 0.3, 0.3), 0, 35),
)


kwargs_text_list_as_tuple_mock = (
    (85.0, 110, 22, "CH3"),
    (185.0, 110, 34, "Y"),
    (253.0, 110, 21, "aMeAla"),
    (385.0, 110, 34, "Q"),
    (485.0, 110, 34, "G"),
    (585.0, 110, 34, "T"),
    (685.0, 110, 34, "F"),
    (785.0, 110, 34, "T"),
    (885.0, 110, 34, "S"),
    (945.0, 165, 34, "D"),
    (885.0, 220, 34, "Y"),
    (785.0, 220, 34, "S"),
    (685.0, 220, 34, "K"),
    (585.0, 220, 34, "Y"),
    (485.0, 220, 34, "L"),
    (385.0, 220, 34, "D"),
    (285.0, 220, 34, "E"),
    (153.0, 220, 19, "Cys(R1)"),
    (125.0, 275, 34, "A"),
    (185.0, 330, 34, "A"),
    (285.0, 330, 34, "K"),
    (385.0, 330, 34, "D"),
    (485.0, 330, 34, "F"),
    (585.0, 330, 34, "V"),
    (653.0, 330, 19, "Cys(R2)"),
    (785.0, 330, 34, "W"),
    (885.0, 330, 34, "L"),
    (945.0, 385, 34, "L"),
    (885.0, 440, 34, "D"),
    (785.0, 440, 34, "H"),
    (685.0, 440, 34, "H"),
    (585.0, 440, 34, "P"),
    (485.0, 440, 34, "S"),
    (385.0, 440, 34, "S"),
    (285.0, 440, 34, "G"),
    (185.0, 440, 34, "Q"),
    (125.0, 495, 34, "P"),
    (185.0, 550, 34, "P"),
    (285.0, 550, 34, "P"),
    (385.0, 550, 34, "S"),
    (478.0, 550, 31, "CN"),
    (585.0, 550, 34, "C"),
    (653.0, 550, 16, "123456789"),
    (755.0, 550, 22, "NH2"),
)


fragments_mock = [
    (["CH3", "Y", "aMeAla", "Q", "G", "T", "F", "T", "S"], 9, "forward", False),
    (["D"], 1, "reverse", True),
    (["Y", "S", "K", "Y", "L", "D", "E", "Cys(R1)"], 8, "reverse", False),
    (["A"], 1, "forward", True),
    (["A", "K", "D", "F", "V", "Cys(R2)", "W", "L"], 8, "forward", False),
    (["L"], 1, "reverse", True),
    (["D", "H", "H", "P", "S", "S", "G", "Q"], 8, "reverse", False),
    (["P"], 1, "forward", True),
    (["P", "P", "S", "CN", "C", "123456789", "NH2"], 7, "forward", False),
]


fragment_kwargs_as_tuple_mock = (
    (70, 100, (0.64453125, 0.828125, 0.64453125), (0.3, 0.3, 0.3), 2, 50),
    (70, 200, (0.58203125, 0.1953125, 0.35546875), (0.3, 0.3, 0.3), 2, 50),
    (70, 300, (0.83984375, 0.2578125, 0.2578125), (0.3, 0.3, 0.3), 2, 50),
    (70, 400, (0.77734375, 0.19921875, 0.19921875), (0.3, 0.3, 0.3), 2, 50),
    (70, 500, (0.390625, 0.61328125, 0.625), (0.3, 0.3, 0.3), 2, 50),
    (70, 600, (0.859375, 0.77734375, 0.3125), (0.3, 0.3, 0.3), 2, 50),
    (70, 700, (0.3984375, 0.2890625, 0.44921875), (0.3, 0.3, 0.3), 2, 50),
    (70, 800, (0.41015625, 0.54296875, 0.41015625), (0.3, 0.3, 0.3), 2, 50),
    (70, 900, (0.109375, 0.25, 0.48828125), (0.3, 0.3, 0.3), 2, 50),
)


fragment_kwargs_txt_as_tuple_mock = (
    (100, 80, 34, "A"),
    (185, 80, 34, "C"),
    (285, 80, 34, "D"),
    (385, 80, 34, "E"),
    (485, 80, 34, "F"),
    (585, 80, 34, "G"),
    (685, 80, 34, "H"),
    (785, 80, 34, "I"),
    (885, 80, 34, "K"),
)

mock_fragment = (["A", "C", "D", "E", "F", "G", "H", "I", "K"], 9, "forward", False)


def as_tuple(kwargs_list, keys=[]):
    """
    Convert a list of dictionaries into a tuple of tuples based on specified keys.
    Args:
        kwargs_list (list): A list of dictionaries to be converted.
        keys (list, optional): A list of keys to extract values from each dictionary. Defaults to an empty list.
    Returns:
        tuple: A tuple of tuples, where each inner tuple contains values from the
         dictionaries corresponding to the specified keys.
    """
    return tuple([tuple([i[key] for key in keys]) for i in kwargs_list])


def test_get_kwargs_from_symbols():
    """
    Test the `get_kwargs_from_symbols` function to ensure it correctly processes
    symbols and returns the expected keyword arguments for drawing.
    The test checks:
    - If the function correctly handles the presence of termini symbols ("N" and "C").
    - If the returned keyword arguments list matches the expected mock data.
    - If the returned text keyword arguments list matches the expected mock data.
    The function `get_kwargs_from_symbols` is expected to return two lists of
     keyword arguments, one for drawing symbols and one for drawing text, which
     are then compared against predefined mock data to validate correctness.
    """
    termini_present = ["N", "C"]
    kwargs_list, kwargs_text_list = get_kwargs_from_symbols(
        copy.deepcopy(problematic_symbols), termini_present=termini_present
    )

    keys = [
        "y",
        "x",
        "rgb_fractions",
        "outline_rgb_fractions",
        "outline_width",
        "radius",
    ]
    kwargs_list_as_tuple = as_tuple(kwargs_list, keys=keys)
    assert kwargs_list_as_tuple == kwargs_list_as_tuple_mock

    txt_keys = ["x", "y", "font_size", "text"]
    kwargs_txt_list_as_tuple = as_tuple(kwargs_text_list, keys=txt_keys)
    assert kwargs_txt_list_as_tuple == kwargs_text_list_as_tuple_mock


def test_generate_kwargs_for_text_in_ellipse_balls():
    """
    Test the function `generate_kwargs_for_text_in_ellipse_balls` to ensure it generates
    the correct keyword arguments for text in ellipse balls.
    This test verifies that the function correctly processes the input parameters and
     produces the expected output dictionary, which is then converted to a tuple for
     comparison with a mock tuple.
    The test checks:
    - The sequence fragment and its length.
    - The direction of the fragment (forward or not).
    - Whether the fragment is a corner.
    - The generated keyword arguments for text in ellipse balls.
    - The conversion of the keyword arguments to a tuple.
    - The equality of the generated tuple with a mock tuple.
    The mock data used in this test includes:
    - `mock_fragment`: A tuple containing the sequence fragment, its length, direction, and corner status.
    - `fragment_kwargs_txt_as_tuple_mock`: The expected tuple of keyword arguments for comparison.
    """
    seq_fragment, fragment_length, fragment_direction, is_corner = mock_fragment
    forward = fragment_direction == "forward"

    kwargs_for_text_in_ellipse_balls = generate_kwargs_for_text_in_ellipse_balls(
        symbols=seq_fragment,
        y=70 + 10,
        forward=forward,
        is_corner=is_corner,
        is_start=True,
        left_margin=100,
        step_x=100,
    )

    kwargs_for_text_in_ellipse_balls_as_tuple = as_tuple(
        kwargs_for_text_in_ellipse_balls, keys=["x", "y", "font_size", "text"]
    )
    assert (
        kwargs_for_text_in_ellipse_balls_as_tuple == fragment_kwargs_txt_as_tuple_mock
    )

    return


def test_generate_kwargs_for_text_in_ellipse_balls2():
    """
    Test the `generate_kwargs_for_ellipse_balls` function to ensure it generates the correct keyword arguments
    for drawing text in ellipse balls.
    This test uses a mock fragment to simulate the input and compares the generated keyword arguments
    with the expected values.
    Steps:
    1. Extract the sequence fragment, fragment length, fragment direction, and corner status from the mock fragment.
    2. Generate the keyword arguments for the ellipse balls using the `generate_kwargs_for_ellipse_balls` function.
    3. Convert the generated keyword arguments to a tuple for comparison.
    4. Assert that the generated keyword arguments match the expected mock values.
    Returns:
        None
    """
    seq_fragment, fragment_length, fragment_direction, is_corner = mock_fragment
    kwargs_in_ellipse_balls = generate_kwargs_for_ellipse_balls(
        symbols=seq_fragment, y=70, forward=True, is_corner=False, is_start=True
    )

    kwargs_in_ellipse_balls_as_tuple = as_tuple(
        kwargs_in_ellipse_balls,
        keys=[
            "y",
            "x",
            "rgb_fractions",
            "outline_rgb_fractions",
            "outline_width",
            "radius",
        ],
    )
    assert kwargs_in_ellipse_balls_as_tuple == fragment_kwargs_as_tuple_mock

    return


def test_get_fragment_kwargs():
    """
    Test the `get_fragment_kwargs` function to ensure it returns the correct
    keyword arguments for drawing a sequence fragment and its corresponding text.
    The test checks:
    - The `fragment_kwargs` dictionary contains the correct values for drawing the fragment.
    - The `fragment_kwargs_txt` dictionary contains the correct values for drawing the fragment's text.
    The function uses mock data for the sequence fragment, its length, direction, and corner status.
    It then converts the returned dictionaries to tuples and compares them with expected mock tuples.
    Assertions:
    - `fragment_kwargs_as_tuple` matches the expected `fragment_kwargs_as_tuple_mock`.
    - `fragment_kwargs_txt_as_tuple` matches the expected `fragment_kwargs_txt_as_tuple_mock`.
    """
    seq_fragment, fragment_length, fragment_direction, is_corner = mock_fragment

    fragment_kwargs, fragment_kwargs_txt = get_fragment_kwargs(
        symbols=seq_fragment,
        y=70,
        is_start=True,
        fragment_direction=fragment_direction,
        is_corner=is_corner,
    )

    fragment_kwargs_as_tuple = as_tuple(
        fragment_kwargs,
        keys=[
            "y",
            "x",
            "rgb_fractions",
            "outline_rgb_fractions",
            "outline_width",
            "radius",
        ],
    )
    fragment_kwargs_txt_as_tuple = as_tuple(
        fragment_kwargs_txt, keys=["x", "y", "font_size", "text"]
    )

    assert fragment_kwargs_as_tuple == fragment_kwargs_as_tuple_mock
    assert fragment_kwargs_txt_as_tuple == fragment_kwargs_txt_as_tuple_mock


def test_get_rev_x_and_font_size():
    """
    Test the `get_rev_x_and_font_size` function with different inputs and font size options.
    The function is tested with the following cases:
    - Single character "A" with variable font size enabled and disabled.
    - String "CH3" with variable font size enabled and disabled.
    The expected results are:
    - For input "A" with variable font size enabled: (15, 34)
    - For input "A" with variable font size disabled: (15, 34)
    - For input "CH3" with variable font size enabled: (30, 28)
    - For input "CH3" with variable font size disabled: (30, 22)
    """
    assert get_rev_x_and_font_size("A", variable_font_size=True) == (15, 34)
    assert get_rev_x_and_font_size("A", variable_font_size=False) == (15, 34)
    assert get_rev_x_and_font_size("CH3", variable_font_size=True) == (30, 28)
    assert get_rev_x_and_font_size("CH3", variable_font_size=False) == (30, 22)


def test_get_start_x():
    """
    Test the `get_start_x` function to ensure it returns the correct starting x-coordinate.
    The test checks the following scenario:
    - left_margin: 100
    - is_corner: False
    - forward: True
    - step_x: 100
    - is_start: True
    Expected result: 100
    """
    assert (
        get_start_x(
            left_margin=100, is_corner=False, forward=True, step_x=100, is_start=True
        )
        == 100
    )


def test_get_N_terminus_params():
    """
    Test the get_N_terminus_params function to ensure it correctly adjusts the
     'x' coordinate of the input parameters dictionary.
    The function is expected to take a dictionary of parameters and return a
     new dictionary with the 'x' coordinate increased by 15 units.
    The input parameters dictionary contains:
        - y: The y-coordinate (int)
        - x: The x-coordinate (int)
        - rgb_fractions: A tuple representing RGB fractions (tuple of floats)
        - outline_rgb_fractions: A tuple representing outline RGB fractions (tuple of floats)
        - outline_width: The width of the outline (int)
        - radius: The radius (int)
    The expected output is a dictionary with the same keys and values, except
     the 'x' coordinate should be increased by 15 units.
    Asserts:
        The function's output matches the expected dictionary with the adjusted 'x' coordinate.
    """
    params = {
        "y": 210,
        "x": 115,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
        "outline_width": 0,
        "radius": 35,
    }
    assert get_N_terminus_params(params) == {
        "y": 210,
        "x": 130,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
        "outline_width": 0,
        "radius": 35,
    }


def test_get_C_terminus_params():
    """
    Test the `get_C_terminus_params` function to ensure it correctly calculates
    the parameters for the C-terminus based on the given `last_params` and
     `previous_params`.
    The test verifies that the function returns the expected dictionary of
     parameters for the C-terminus, including position (`x`, `y`), color
     (`rgb_fractions`, `outline_rgb_fractions`), outline width (`outline_width`),
     and radius (`radius`).
    The expected output is compared against the actual output using an assertion.
    """
    last_params = {
        "y": 430,
        "x": 185.0,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
        "outline_width": 0,
        "radius": 35,
    }

    previous_params = {
        "y": 375,
        "x": 140.0,
        "rgb_fractions": (0.37109375, 0.58203125, 0.51171875),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
        "outline_width": 2,
        "radius": 50,
    }
    c_terminus_params = get_C_terminus_params(
        params=last_params, previous_params=previous_params
    )

    assert c_terminus_params == {
        "y": 430,
        "x": 170.0,
        "rgb_fractions": (0.3, 0.3, 0.3),
        "outline_rgb_fractions": (0.3, 0.3, 0.3),
        "outline_width": 0,
        "radius": 35,
    }


def test_draw_sequence():
    """
    Test the draw_symbols function with various sequences of symbols.
    This test function creates multiple lists of symbols representing sequences
    and calls the draw_symbols function to draw them with specified dimensions.
    It tests the function with different sequences to ensure it handles various
    inputs correctly.
    The sequences include:
    - A sequence with cysteine residues and other amino acids.
    - A sequence with modified amino acids and other residues.
    - A sequence with additional symbols and a longer length.
    The function also tests drawing the same sequence twice, once with a deep copy
    and once with an output file specified.
    Returns:
        None
    """
    """
    test drawing sequence from symbols
    """
    symbols = [
        "Cys(R1)",
        "S",
        "Cys(R2)",
        "S",
        "D",
        "E",
        "K",
        "S",
        "D",
        "E",
        "K",
        "M",
        "N",
        "P",
        "Q",
    ] + [
        "Cys(R1)",
        "S",
        "Cys(R2)",
        "S",
        "D",
        "E",
        "K",
        "S",
        "D",
        "E",
        "K",
        "M",
        "N",
        "P",
        "Q",
        "R",
    ]
    draw_symbols(symbols, width=1024, height=1024)

    symbols = (
        [
            "CH3",
            "Y",
            "aMeAla",
            "Q",
            "G",
            "T",
            "F",
            "T",
            "S",
            "D",
            "Y",
            "S",
            "K",
            "Y",
            "L",
        ]
        + [
            "D",
            "E",
            "Cys(R1)",
            "A",
            "A",
            "K",
            "D",
            "F",
            "V",
            "Cys(R2)",
            "W",
            "L",
            "L",
            "D",
            "H",
            "H",
        ]
        + ["P", "S", "S", "G", "Q", "P", "P", "P", "S", "NH2"]
    )
    draw_symbols(symbols, width=1024, height=1024)

    symbols = (
        [
            "CH3",
            "Y",
            "aMeAla",
            "Q",
            "G",
            "T",
            "F",
            "T",
            "S",
            "D",
            "Y",
            "S",
            "K",
            "Y",
            "L",
        ]
        + [
            "D",
            "E",
            "Cys(R1)",
            "A",
            "A",
            "K",
            "D",
            "F",
            "V",
            "Cys(R2)",
            "W",
            "L",
            "L",
            "D",
            "H",
            "H",
        ]
        + ["P", "S", "S", "G", "Q", "P", "P", "P", "S", "CN", "C", "123456789", "NH2"]
    )
    draw_symbols(copy.deepcopy(symbols), width=1024, height=1024)
    draw_symbols(
        symbols, width=1024, height=1024, out="../image_file_name_can_err_2.png"
    )
    return
