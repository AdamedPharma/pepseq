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


*****************************************************
Reading Sequence and Modification from SMILES string
*****************************************************


To calculate sequence(s) in Pepseq Format and  list of modification SMILES from
Modified Peptide structure given in SMILES you can use the
``pepseq.commands.read_smiles()`` function:

.. autofunction:: commands.read_smiles

.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.read-smiles('stuff_in.smi', db_path='augmented_db.json', out=path)

************************************
Augment Monomer Database
************************************

To augment monomers database with new monomers and modifications you can use
the ``pepseq.commands.augment_db_json_command()``

.. autofunction:: commands.augment_db_json_command

.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.augment_db_json_command(sdf_path='sdf_file.sdf',
                                out='db_json_augmented.json'
                                )


********************************************
Calculate JSON representing Modified Peptide
********************************************

to calculate JSON representing Modified Peptide through peptide sequence in
Pepseq format and modification list in SMILES format
you can use the ``pepseq.commands.calculate_json_from()`` function:

.. autofunction:: commands.calculate_json_from

.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.calculate_json_from
        ('{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF',
            ['[1*]CNCC[2*]' --mod-smiles '[3*]CNCCSP']
