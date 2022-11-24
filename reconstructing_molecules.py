import rdkit

"""

Jesli to jest cos co chcemy jeszcze potem w przyszlosci deponowac gdzies jako open
source to faktycznie docelowo i tak ma miec osobne repozytorium,  wiec mozemy
mu dac osobne repozytorium juz teraz ...
Kod jako API to z jednej strony pewnie to powinno dzialac w linii polecen, ale tez
powinno byc do zaimportowania jakos do LabGears  jako funkcja ...

Jesli chodzi o linie polecen to widze to tak np. ze Python
 pepseq_commands.py generate_smiles 'IAMSEQUENCE'  '[2*]CCCNNCI'  '[5]CCCNNC[7*]I'
     i zwraca 'CCCIAMSMILESOFWHOLEPEPTIDE'

no i ten singledispatch moze byc jakos wewnatrz, zeby bardziej elegancko parametry
 ogarnac do command line: CLI w typer(ze)

no i tam moze byc jakaa klasa o nazwie np. ModifiedPeptide ktora bedzie _init_ -owana
jakas metoda dekorowana singledispatch-em
printy beda chyba musialy miec dwa formaty ...

Bo chyba arbitralnych numerow - tak zeby pasowaly do nr residue sie, nie da? ...
  oczywiscie ze sie da i to jest dobry pomysl

"""
    


class UseCaseCombineCombo(object):

    def __init__(self, combo=None):
        self.combo = combo
        return

    @property   
    def edcombo(self):
        if self.__dict__.get('_edcombo') is None:
            self._edcombo = rdkit.Chem.EditableMol(self.combo)
        return self._edcombo

    def radicals(self, radical_name):
        rads = []
        for atom in self.combo.GetAtoms():
            if atom.GetIsotope() == radical_name:
                rads.append( atom )
        return rads

    def execute(self):
        """

        simplest case:

            N-S-R1 is to be connected with C-O-R1

            S is connected with O

            Resulting Molecule is N-S-O-C

            Combo containing:

            N <- 1
            S <- 2
            R1(ns) <- 3
            C <- 4
            O <- 5
            R1(co) <- 6

            S2 <-(is_connected_to)-> O5

            S2 <-(is_neighbor_of)-> R1(ns)3

            O5 <-(is_neighbor_of)-> R1(co)6

        more advanced case:

        R2-N-S-R1 is to be connected with R2-C-O-R1

        S is connected with O (R1)
        N is connected with C (R2)

        Resulting Molecule is:

         N-S
         | |
         C-O
        
        """

        radical_names = self.radical_names

        for radical_name in radical_names:
            self.bond_radical( radical_name )
            self.remove_radical( radical_name )
            mol = self.edcombo.GetMol()
            self.combo = mol

        return mol
    
    def remove_radical(self, radical_name):
        radicals = self.radicals(radical_name)
        radical_ids = [ radical.GetIdx() for radical in radicals ]
        for radical_id in reversed( sorted(  radical_ids ) ):
            self.edcombo.RemoveAtom( radical_id )
        return

    def bond_radical(self, radical_name):
        radicals = self.radicals(radical_name)

        num_radicals = len( radicals )
        
        for i in range( num_radicals ):
            radical_i = radicals[ i ]
            neighbors_i = radical_i.GetNeighbors()

            for j in range( i+1, num_radicals):
                radical_j = radicals[ j ]
                neighbors_j = radical_j.GetNeighbors()

                #we join neighbors instead of radicals
                for neighbor_i in neighbors_i:
                    for neighbor_j in neighbors_j:
                        self.edcombo.AddBond( neighbor_i.GetIdx(), neighbor_j.GetIdx() )
        return

    def get_radical_names(self):
        """
        radical_names = set( [1, 2 ] )
        """

        radical_names = set( [ ] )

        for atom in self.combo.GetAtoms():
            isotope = atom.GetIsotope()
            if isotope != 0:
                radical_names.add( isotope )
        return radical_names
    
    @property
    def radical_names(self):
        if self.__dict__.get('_radical_names') is None:
            self._radical_names = self.get_radical_names()
        return self._radical_names

