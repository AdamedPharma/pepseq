.. _commands:

*******************************
draw
*******************************

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

************************************
draw_pepseq
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.draw.draw_pepseq