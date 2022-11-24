# pepseq
Cheminformtic library to work flawlessly with modified peptide structures

Aim: A text format to conveniently handle complicated modified peptide structures as standard formats like SMILES or Fasta sequence are either to simple or to complicated.

New format is a combination of:

Easily readable 1aa code like : 
  o	KYLDERAAQDFVQW

3aa code that allows to handle non-standard amino acids
  o	Ala-Aib-Gly-Lys
Common chemical names to define typical chemical substitution of amino acids, like: 
o	Thr(tBu), Ser(tBu), Phe(4-Cl)

SMILES codes to describe unspecified modifications:
KYLDCRAAQDFVQW<[5*]CCCC(N)=O>
