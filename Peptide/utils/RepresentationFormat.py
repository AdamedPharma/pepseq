class RepresentationFormat(object):
    rules = {"mask_mod_aa": None, "show_termini": True}

    def __init__(self, name: str = None):
        self.name = name
        return

    def json():
        """
        e. g.
        {
            'num_letters': 3,
            'separator': '-',
            'Aib':'J'
            }

        """
        return

    def write():
        return  # here? not in parser?

    def reads_sequence_txt():
        return  # here? not in parser?


class AlignmentRepresentation(RepresentationFormat):
    """ """

    rules = {"mask_mod_aa": "X", "show_termini": False}

    def __init__(self, name: str = "AlignmentRepresentation"):
        self.name = name
        return

    def x(self):
        return


class DallasMonomersRepresentation(RepresentationFormat):
    """
        Rule #1: use of D- and DL- prefixes
            Widely understood and accepted stereo configuration prefixes 'D-', 'L-' (default) and 'DL-' should be used.
            example: D-Arg; DL-Ala; D-Ser

        Rule #2: retain widely used 3LC's.
        Universally accepted code/definitions should be used:
            Abu (2-aminobutyric acid), Aib (2-aminoisobutyric acid), aIle (allo-isoleucine),
            aThe (allo-threonine), bAla (beta-alanine), Cha (beta-cyclohexylalanine),
            Chg (alpha-cyclohexylglycine), Cit (Citrulline), Dab (2,3-diaminobutyric acid),
            Dap (2,3-diaminopropionic acid), hArg (homoarginine), Hcy (homocysteine),
            Hse (homoserine), hPhe (homophenylalanine), Nle (norleucine), Nva (norvaline),
            Orn (ornithine), Pen (penicillamine), Phg (phenylglycine), Sar (sarcosine),
            Sec (selenocysteine), 2Thi (thien-2-ylalanine), 3Thi (thien-3-ylalanine),
            xiIle (xi-isoleucine), xiThr (xi-threonine), and many more.

        Rule #3: modifying prefixes
            'N(Me)': N-modified variants. 'N(Me)Gly' = 'Sar'.
            'O': depsi-peptides that connect via ester linkages, e.g. 'Arg-OAla-Ser'.
            'aMe': alpha-methyl variants (aMeAla = Aib)

        Rule #4: line formulae substs
            Substitutions are indicated by line formulae
            Ph - phenyl; CF3 - trifluoromethyl; Me - methyl; Tos - tosyl; Ac - acetyl;
            Trt - trityl; Cl - chloro; Tf - triflyl; CN - cyano; SO3H - sulfonic acid;
            Bn - benzyl; OPfp - perfluorophenoxy; NH2 - amino; iPr - isopropyl;
            NO2 - nitro; For - formyl; tBu - tert-butyl ...

        Rule #5: default substitution locants
            Default substitution locants: Abu - CG; Ala - CB; Arg - NH2; Cys - SG;
            Dab - ND; Dap - NG; Gly - CA; Lys - NZ; Ser - OG; Thr - OG1; Tyr - OH

        Rule #6: specified substitution locants
            Amino acids may be substituted at particular locant using the following syntax:
                Phe(4-Cl); Phg(4-CH2NH2); Tyr(2,6-diCl-Bn);
            Substitutions may occur at more than one locant:
                His(2,5-diI);
            Tyr(2,6-diF) <- dajemy sobie spokoj z podstawieniami pierscienia


        Rule #7: implicit leaving groups
            7.1 hydrogen for most Amino Acids (e.g. Ala(Cl) )
            7.2 OH for carboxylic Amino Acids (e.g. Asp(NH2) ).
                Requires "oxy" to form esters (Asp(OMe))

            7.3. Cysteine (and penicillamine) substitute on the sulfur (isn't that already in Rule #5):
                Cys(tBu); Cys(StBu)

        Examples:

        1. H~CPY~OH, H-Cys-Pro-Trp-His-Leu-Leu-Pro-Phe-Cys-OH, CHEMBL501567, CC(C)C[C@H](NC(=O)[C@H](Cc1cnc[nH]1)NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](N)CS)C(=O)N[C@@H](CC(C)C)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CS)C(=O)O
        2. H-Tyr-Pro-Phe-Phe-OtBu, CHEMBL500195

        3. cyclo[Ala-Tyr-Val-Orn-Leu-D-Phe-Pro-Phe-D-Phe-Asn], CHEMBL438006

        4. H-Nle(Et)-Tyr-Pro-Trp-Phe-NH2, CHEMBL500704

        5. H-{SMILES tego Nle(Et)}-Tyr-Pro-Trp-Phe-NH2, CHEMBL500704

        6. H-DL-hPhe-Val-Met-Tyr(PO3H2)-Asn-Leu-Gly-Glu-OH, CHEMBL439086 (najbardziej straightforward) (tyr najlepszzy)

        7. cyclo[Phe-D-Trp-Tyr(Me)-D-Pro], CHEMBL507127

        7.1 cyclo[Phe-D-Trp-Tyr('C')-D-Pro], CHEMBL507127


        8.  H-D-Pyr-D-Leu-pyrrolidide, CHEMBL1181307

        9. Ac-DL-Phe-aThr-Leu-Asp-Ala-Asp-DL-Phe(4-Cl)-OH, CHEMBL1791047

        10. H-D-Cys(1)-D-Asp-Gly-Tyr(3-NO2)-Gly-Hyp-Asp-D-Cys(1)-NH2, CHEMBL583516

        11. Boc-Tyr-Tyr(3-Br)-OMe, CHEMBL1976073

        12. CC(C)[C@H]1C(=O)N[C@H](C(=O)N[C@H]( C(=O)N2CCC[C@H]2C(=O)N[C@H](C(=O)N [C@H](C(=O)N1)CC3=CC=C(C=C3)OP(=O)( O)O)CSSC[C@@H](C(=O)O)NC(=O)C)C(C) C)CC(=O)N, cyclo[Asn-Val-Pro-Cys(1)-Tyr(PO3H2)-Val].
    Ac-Cys(1)-OH
        testy z baza nasza

        13. test na oznaczenie zmodyfikowanego w kolku (oznaczyz tylko ze aminokwasy jest zmodyfikowany: i oznaczyc zrodlowy)
        problemy z kolejnoscia w cyklicznych aminokwasach


    """


class DallasMonomersRepresentation2(RepresentationFormat):
    """
    Rule #1: use of D- and DL- prefixes
        Widely understood and accepted stereo configuration prefixes 'D-', 'L-' (default) and 'DL-' should be used.
        example: D-Arg; DL-Ala; D-Ser
    """

    """
    Rule #2: retain widely used 3LC's.
    Universally accepted code/definitions should be used: 
        Abu (2-aminobutyric acid), Aib (2-aminoisobutyric acid), aIle (allo-isoleucine),
        aThe (allo-threonine), bAla (beta-alanine), Cha (beta-cyclohexylalanine),
        Chg (alpha-cyclohexylglycine), Cit (Citrulline), Dab (2,3-diaminobutyric acid),
        Dap (2,3-diaminopropionic acid), hArg (homoarginine), Hcy (homocysteine),
        Hse (homoserine), hPhe (homophenylalanine), Nle (norleucine), Nva (norvaline),
        Orn (ornithine), Pen (penicillamine), Phg (phenylglycine), Sar (sarcosine),
        Sec (selenocysteine), 2Thi (thien-2-ylalanine), 3Thi (thien-3-ylalanine),
        xiIle (xi-isoleucine), xiThr (xi-threonine), and many more.
    """

    """
    Rule #3: modifying prefixes 
        'N(Me)': N-modified variants. 'N(Me)Gly' = 'Sar'.
        'O': depsi-peptides that connect via ester linkages, e.g. 'Arg-OAla-Ser'. 
        'aMe': alpha-methyl variants (aMeAla = Aib)
    """

    """
    Rule #4: line formulae substs
        Substitutions are indicated by line formulae 
        Ph - phenyl; CF3 - trifluoromethyl; Me - methyl; Tos - tosyl; Ac - acetyl;
        Trt - trityl; Cl - chloro; Tf - triflyl; CN - cyano; SO3H - sulfonic acid;
        Bn - benzyl; OPfp - perfluorophenoxy; NH2 - amino; iPr - isopropyl; 
        NO2 - nitro; For - formyl; tBu - tert-butyl ...
    """

    """
    Rule #5: default substitution locants
        Default substitution locants: Abu - CG; Ala - CB; Arg - NH2; Cys - SG; 
        Dab - ND; Dap - NG; Gly - CA; Lys - NZ; Ser - OG; Thr - OG1; Tyr - OH
    """

    """
    Rule #6: specified substitution locants
        Amino acids may be substituted at particular locant using the following syntax: 
            Phe(4-Cl); Phg(4-CH2NH2); Tyr(2,6-diCl-Bn); 
        Substitutions may occur at more than one locant:
            His(2,5-diI); 
        Tyr(2,6-diF) <- dajemy sobie spokoj z podstawieniami pierscienia
    """

    """
    Rule #7: implicit leaving groups
        7.1 hydrogen for most Amino Acids (e.g. Ala(Cl) )
        7.2 OH for carboxylic Amino Acids (e.g. Asp(NH2) ).
            Requires "oxy" to form esters (Asp(OMe))
     
        7.3. Cysteine (and penicillamine) substitute on the sulfur (isn't that already in Rule #5):
            Cys(tBu); Cys(StBu)

    """

    """
    Examples:
            
    1. H~CPY~OH, H-Cys-Pro-Trp-His-Leu-Leu-Pro-Phe-Cys-OH, CHEMBL501567, CC(C)C[C@H](NC(=O)[C@H](Cc1cnc[nH]1)NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](N)CS)C(=O)N[C@@H](CC(C)C)C(=O)N1CCC[C@H]1C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CS)C(=O)O
    2. H-Tyr-Pro-Phe-Phe-OtBu, CHEMBL500195 

    3. cyclo[Ala-Tyr-Val-Orn-Leu-D-Phe-Pro-Phe-D-Phe-Asn], CHEMBL438006 

    4. H-Nle(Et)-Tyr-Pro-Trp-Phe-NH2, CHEMBL500704 

    5. H-{SMILES tego Nle(Et)}-Tyr-Pro-Trp-Phe-NH2, CHEMBL500704 

    6. H-DL-hPhe-Val-Met-Tyr(PO3H2)-Asn-Leu-Gly-Glu-OH, CHEMBL439086 (najbardziej straightforward) (tyr najlepszzy)

    7. cyclo[Phe-D-Trp-Tyr(Me)-D-Pro], CHEMBL507127 

    7.1 cyclo[Phe-D-Trp-Tyr('C')-D-Pro], CHEMBL507127 


    8.  H-D-Pyr-D-Leu-pyrrolidide, CHEMBL1181307 

    9. Ac-DL-Phe-aThr-Leu-Asp-Ala-Asp-DL-Phe(4-Cl)-OH, CHEMBL1791047 

    10. H-D-Cys(1)-D-Asp-Gly-Tyr(3-NO2)-Gly-Hyp-Asp-D-Cys(1)-NH2, CHEMBL583516 

    11. Boc-Tyr-Tyr(3-Br)-OMe, CHEMBL1976073

    12. CC(C)[C@H]1C(=O)N[C@H](C(=O)N[C@H]( C(=O)N2CCC[C@H]2C(=O)N[C@H](C(=O)N [C@H](C(=O)N1)CC3=CC=C(C=C3)OP(=O)( O)O)CSSC[C@@H](C(=O)O)NC(=O)C)C(C) C)CC(=O)N, cyclo[Asn-Val-Pro-Cys(1)-Tyr(PO3H2)-Val]. 
Ac-Cys(1)-OH
    testy z baza nasza

    13. test na oznaczenie zmodyfikowanego w kolku (oznaczyz tylko ze aminokwasy jest zmodyfikowany: i oznaczyc zrodlowy)
    problemy z kolejnoscia w cyklicznych aminokwasach

    
    """
