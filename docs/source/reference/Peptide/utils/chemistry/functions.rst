.. _commands:

*******************************
functions
*******************************

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

************************************
get_connecting_bonds
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.Peptide.utils.chemistry.functions.get_connecting_bonds
