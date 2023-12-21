.. _commands:

functions
=========

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

validate
--------

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.functions.validate()`` function:


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.functions.validate('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.functions.validate

.. autofunction:: pepseq.functions.calculate_pepseq_and_mods

.. autofunction:: pepseq.functions.calculate
