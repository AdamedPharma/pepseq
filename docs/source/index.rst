.. _pepseq_docs_mainpage:

####################
Pepseq documentation
####################



.. toctree::
   :maxdepth: 1
   :hidden:

   User Guide <user/index>
   API reference <api>
   Code reference <reference/index>


**Version**: |version|

**Download documentation**:
`Historical versions of documentation <https://readthedocs.org/pepseq/doc/>`_
   
**Useful links**:
`Source Repository <https://github.com/pepseq>`_ |
`Issue Tracker <https://github.com/pepseq/issues>`_ |

.. grid:: 2


    .. grid-item-card::
        :img-top: ../source/_static/index-images/user_guide.svg


        Code reference
        ^^^^^^^^^^^^^

        The reference guide contains a detailed description of the functions,
        modules, and objects included in NumPy. The reference describes how the
        methods work and which parameters can be used. It assumes that you have an
        understanding of the key concepts.

        +++

        .. button-ref:: reference/index
            :expand:
            :color: secondary
            :click-parent:

            To the reference guide


**Pepseq** is a Python library for chemoinformaticians that
creates a new format for representing modified peptide called Pepseq.
Format can be created from input arguments or from parsed *SMILES* code.

To build this documentation cd to main project folder and type:
sphinx-build -M html docs/source/ docs/build/

.. note::

   This project is under active development.




