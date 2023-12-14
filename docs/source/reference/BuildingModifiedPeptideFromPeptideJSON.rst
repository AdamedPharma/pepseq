.. _Building_Modified_Peptide_From_PeptideJSON:

*******************************
Building Modified Peptide From PeptideJSON
*******************************

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

************************************
Add Internal Bond
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.add_internal_bond

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.add_disulfide_bond

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_attachment_points

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.find_atom

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.add_staple

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_peptide_json_from_sequence
   
.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_residue_symbols_from_sequence

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_molecule_from_sequence

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_smiles_from_sequence

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_molecule_from_json

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_smiles_from_peptide_json
