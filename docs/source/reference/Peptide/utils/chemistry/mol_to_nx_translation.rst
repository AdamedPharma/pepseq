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


.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_edge_tuple

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.mol_to_nx
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.nx_to_mol

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_chiral_tag_int

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_hybridization_int

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_node_tuple
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.nx_to_json

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_mol_json
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.mol_json_to_nx
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.mol_json_to_mol
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.draw_peptide_json
