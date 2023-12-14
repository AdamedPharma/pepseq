.. _commands:

*******************************
Command Line Interface Commands
*******************************

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

************************************
Creating SMILES from Pepseq Format
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.Drawer.get_start_x

.. autofunction:: pepseq.Drawer.get_rev_x_and_font_size

.. autofunction:: pepseq.Drawer.generate_kwargs_for_text_in_ellipse_balls

.. autofunction:: pepseq.Drawer.generate_kwargs_for_ellipse_balls

.. autofunction:: pepseq.Drawer.get_is_corner

.. autofunction:: pepseq.Drawer.get_direction

.. autofunction:: pepseq.Drawer.get_fragment_length

.. autofunction:: pepseq.Drawer.schema_layout_generator_from_symbols

.. autofunction:: pepseq.Drawer.get_fragment_kwargs

.. autofunction:: pepseq.Drawer.get_N_terminus_params

.. autofunction:: pepseq.Drawer.get_C_terminus_params

.. autofunction:: pepseq.Drawer.get_kwargs_from_symbols

.. autofunction:: pepseq.Drawer.draw_ellipse_ball

.. autofunction:: pepseq.Drawer.draw_ellipse_balls

.. autofunction:: pepseq.Drawer.draw_text_in_ellipse

.. autofunction:: pepseq.Drawer.draw_text_in_ellipse_balls

.. autofunction:: pepseq.Drawer.get_png_string_from_surface

.. autofunction:: pepseq.Drawer.draw_symbols
