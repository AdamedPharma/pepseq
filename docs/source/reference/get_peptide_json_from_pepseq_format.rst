.. _commands:

get_peptide_json_from_pepseq_format
===================================

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

decompose_symbol
----------------

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.get_peptide_json_from_pepseq_format.get_attachment_point_json()`` function:


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'



.. autofunction:: pepseq.Peptide.utils.pure_parsing_functions.decompose_symbol

.. autofunction:: pepseq.Peptide.utils.pure_parsing_functions.get_attachment_points_on_sequence_json

.. autofunction:: pepseq.Peptide.utils.pure_parsing_functions.get_base_seq

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_single_modification_json

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_ext_mod_json
   
.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_pep_json
