.. _commands:

Parser
^^^^^^

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

output_modified_residue
""""""""""""""""""""""""

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.Peptide.utils.Parser.output_modified_residue
   
.. autofunction:: pepseq.Peptide.utils.Parser.parentheses_locs_list
   
.. autofunction:: pepseq.Peptide.utils.Parser.find_parentheses

.. autofunction:: pepseq.Peptide.utils.Parser.parse_canonical

.. autofunction:: pepseq.Peptide.utils.Parser.parse_canonical2
   
.. autofunction:: pepseq.Peptide.utils.Parser.find_termini
   
.. autofunction:: pepseq.Peptide.utils.Parser.get_canonical
