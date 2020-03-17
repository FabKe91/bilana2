'''
    Contains all atom name definitions for common lipid molecules in the
    CHARMM36 force field

    Constants are:
        - GLYCEROLATOMS         List of all atoms in the glycerol backbone
        - HEAD_ATOMS_OF         Dict of all atoms in head group of lipid type
        - TAILCARBONS           Dict of all carbon atoms of tail type
        - TAILHYDROGENS         Dict of all hydrogen atoms of tail type
        - TAIL_ATOMS_OF         Combining tail carbons and tail hydrogens
        - SCD_TAIL_ATOMS_OF     List of the sequence of atom names for SCD calculation
        - CENTRAL_ATOM_OF       Returns the name of the reference atom of respective molecule

        - TAILS                 Lists all known tail types defined in this module
        - HEADS                 Lists all known head types defined in this module
        - STEROLS               Lists all known sterol molecules in this module
        - PROTEIN_SEQUENCES     Lists all known protein sequence names
        - PROTEIN_RESIDUES      Lists all known protein residue types #### In future all residue types should be included
        - SOLVENTS              Lists the names of all solvent types in the respective force field

'''

# =====================================================
# Implemented molecules, chain types, protein sequences
# Add all new definitions here so in the tests you can
# be sure everything was properly added
# =====================================================

TAILS             = ['DP', 'DM', 'DS', 'DO', 'DY', 'DU', 'PO', 'PL']
HEADS             = ['PC', 'PE', 'PS', 'PI', 'PA']
STEROLS           = ['CHL1', 'CHIM', 'CH0M', 'ERG', 'ch1m']
PROTEIN_SEQUENCES = ['WSC1', ]
PROTEIN_RESIDUES  = ['VAL', 'GLY', 'ALA', 'ILE', 'LEU', 'CYS', 'ARG', 'HSD']
IONS              = ["CL", "POT", "NA" ]
WATER             = ["TIP3", "SOL",]
SOLVENTS          = WATER + IONS


# =================================================
# Glycerolpart, also carbonylpart of FA is included
# =================================================

GLYCEROLATOMS = ['C1', 'O11', 'C2', 'O21', 'C21', 'O22', 'C3', 'O31', 'C31', 'O32', 'HA', 'HB', 'HY', 'HX', 'HS', ]

# =================================================
# Head atom definition
# =================================================

HEAD_ATOMS_OF = {
    'PC':['P', 'O12', 'O13', 'O14', 'N', 'C11', 'C12', 'C13', 'C14', 'C15',
          'H11A', 'H11B', 'H12A', 'H12B', 'H13A', 'H13B', 'H13C', 'H14A', 'H14B', 'H14C', 'H15A', 'H15B', 'H15C'],
    'PA':['P', 'O12', 'O13', 'O14', 'H12'],
    'PE':['P', 'O12', 'O13', 'O14', 'C11', 'C12', 'N', 'HN1', 'HN2', 'HN3', 'H11A', 'H11B', 'H12A', 'H12B',],
    'PI':['P', 'O12', 'O13', 'O14', 'C11', 'C12', 'O2', 'C13', 'O3', 'C14', 'O4', 'C15', 'O5', 'C16', 'O6',
          'H1', 'H2', 'HO2', 'H3', 'HO3', 'H4', 'HO4', 'H5', 'HO5', 'H6', 'HO6'],
    'PS':['P', 'O12', 'O13', 'O14', 'C11', 'C12', 'C13', 'O13A', 'O13B', 'N', 'H11A', 'H11B', 'H12A', 'HN1', 'HN2', 'HN3'],
    #### The entries for sterols should be complete *with* hydrogens and [head atoms]+[tail atoms] ###
    #### should be a list of all atoms in a sterol molecule ###
    'CHL1':['O3', 'C1', 'C2', 'C3', 'C4', 'C5', 'C10'],
    'CHIM':['C20', "C1", "C2", "N1", "N2"],
    'CH0M':['C20', "C1", "C2", "N1", "N2"],
    'ERG':['O3', 'C1', 'C2', 'C3', 'C4', 'C5', 'C10'],
    }

# =================================================
# Tail definitions
# =================================================

TAILCARBONS = {
    'DP':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216'],                  #16:0
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:0
    'DM':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214',],                                 #14:0
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314',]],                                #14:0
    'DS':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:0
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:0
    'DU':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:2
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:2
    'DY':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216'],                  #16:1
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:1
    'PL':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:2
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:0
    'PO':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:1
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']],                 #16:0
    'DO':[['C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],  #18:1
          ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316', 'C317', 'C318']], #18:1
    }

TAILHYDROGENS = {
    'DP':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S', 'H7R', 'H7S', 'H8R', 'H8S', 'H9R', 'H9S', 'H10R', 'H10S', 'H11R', 'H11S', 'H12R', 'H12S', 'H13R', 'H13S', 'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],
          ['H2X', 'H2Y','H3X', 'H3Y', 'H4X', 'H4Y', 'H5X', 'H5Y', 'H6X', 'H6Y', 'H7X', 'H7Y', 'H8X', 'H8Y', 'H9X', 'H9Y', 'H10X', 'H10Y', 'H11X', 'H11Y', 'H12X', 'H12Y', 'H13X', 'H13Y', 'H14X', 'H14Y', 'H15X', 'H15Y','H16X', 'H16Y', 'H16Z']],
    'DM':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                'H14R', 'H14S', 'H14T'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y', 'H14Z']],
    'DS':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H9S', 'H10R', 'H10S','H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17R', 'H17S', 'H118R', 'H18S', 'H18T'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y','H15X', 'H15Y','H16X', 'H16Y','H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
    'DU':[['H2R', 'H2S', 'H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S', 'H7R', 'H7S','H8R', 'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H13R', 'H14R', 'H14S', 'H15R', 'H15S', 'H16R', 'H16S', 'H17R', 'H17S', 'H18R', 'H18S', 'H18T'],
          ['H2X', 'H2Y', 'H3X', 'H3Y', 'H4X', 'H4Y', 'H5X', 'H5Y', 'H6X', 'H6Y', 'H7X', 'H7Y','H8X', 'H8Y', 'H9X', 'H10X', 'H11X', 'H11Y', 'H12X', 'H13X', 'H14X', 'H14Y', 'H15X', 'H15Y', 'H16X', 'H16Y', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
    'DY':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H16T'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H10X', 'H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
    'PL':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H13R',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H9Y','H10X', 'H10Y','H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
    'PO':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H91', 'H101', 'H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17R', 'H17S', 'H118R', 'H18S', 'H18T'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H9Y','H10X', 'H10Y', 'H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y','H15X', 'H15Y','H16X', 'H16Y','H16Z']],
    'DO':[['H2R', 'H2S','H3R', 'H3S', 'H4R', 'H4S', 'H5R', 'H5S', 'H6R', 'H6S','H7R','H7S','H8R',
                'H8S', 'H9R', 'H10R', 'H11R', 'H11S', 'H12R', 'H12S','H13R', 'H13S',
                'H14R', 'H14S', 'H15R', 'H15S','H16R', 'H16S', 'H17R', 'H17S', 'H118R', 'H18S', 'H18T'],
             ['H2X', 'H2Y','H3X', 'H3Y','H4X', 'H4Y','H5X', 'H5Y', 'H6X', 'H6Y','H7X', 'H7Y','H8X',
                'H8Y','H9X', 'H10X', 'H11X', 'H11Y','H12X', 'H12Y', 'H13X', 'H13Y','H14X',
                'H14Y','H15X', 'H15Y','H16X', 'H16Y','H17X', 'H17Y', 'H18X', 'H18Y', 'H18Z']],
    }

TAIL_ATOMS_OF = {
    'DP':[TAILCARBONS['DP'][0], TAILHYDROGENS['DP'][0], TAILCARBONS['DP'][1], TAILHYDROGENS['DP'][1],],
    'DM':[TAILCARBONS['DM'][0], TAILHYDROGENS['DM'][0], TAILCARBONS['DM'][1], TAILHYDROGENS['DM'][1],],
    'DS':[TAILCARBONS['DS'][0], TAILHYDROGENS['DS'][0], TAILCARBONS['DS'][1], TAILHYDROGENS['DS'][1],],
    'DO':[TAILCARBONS['DO'][0], TAILHYDROGENS['DO'][0], TAILCARBONS['DO'][1], TAILHYDROGENS['DO'][1],],
    'DY':[TAILCARBONS['DY'][0], TAILHYDROGENS['DY'][0], TAILCARBONS['DY'][1], TAILHYDROGENS['DY'][1],],
    'DU':[TAILCARBONS['DU'][0], TAILHYDROGENS['DU'][0], TAILCARBONS['DU'][1], TAILHYDROGENS['DU'][1],],
    'PO':[TAILCARBONS['PO'][0], TAILHYDROGENS['PO'][0], TAILCARBONS['PO'][1], TAILHYDROGENS['PO'][1],],
    'PL':[TAILCARBONS['PL'][0], TAILHYDROGENS['PL'][0], TAILCARBONS['PL'][1], TAILHYDROGENS['PL'][1],],

    #### The entries for sterols should be complete *with* hydrogens and [head atoms]+[tail atoms] ###
    #### should be a list of all atoms in a sterol molecule ###
    'CHL1':[['C13', 'C14', 'C15', 'C16', 'C17', 'C20', 'C22', 'C23', 'C24', 'C25']],
    'ERG':[['C13', 'C14', 'C15', 'C16', 'C17', 'C20', 'C22', 'C23', 'C24', 'C25']],
    'CHIM':[['C12', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30']],
    'CH0M':[['C12', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30']],
    }

# ====================================================
# In the following special groups of atoms are defined
# ====================================================

SCD_TAIL_ATOMS_OF = {
    'DP':[TAILCARBONS['DP'][0][::2], TAILCARBONS['DP'][1][::2]],
    'DU':[['C22', 'C24', 'C26', 'C28', 'C211', 'C214', 'C216', 'C218'], # Double bonds between 9-10, 12-12
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C314', 'C316', 'C318']],
    'DM':[TAILCARBONS['DM'][0][::2], TAILCARBONS['DM'][1][::2]],
    'DS':[TAILCARBONS['DS'][0][::2], TAILCARBONS['DS'][1][::2]],
    'DY':[['C22', 'C24', 'C26', 'C28', 'C211', 'C213', 'C215', ],   # Double bonds between 9-10
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C313', 'C315',]],
    'PL':[['C22', 'C24', 'C26', 'C28', 'C211', 'C214', 'C216', 'C218'], TAILCARBONS['PL'][1][::2]],
    'PO':[['C22', 'C24', 'C26', 'C28', 'C211', 'C213', 'C215', 'C217'], TAILCARBONS['PO'][1][::2]],
    'DO':[['C22', 'C24', 'C26', 'C28', 'C211', 'C213', 'C215', 'C217'],  # Double bonds between 9-10
            ['C32', 'C34', 'C36', 'C38', 'C311', 'C313', 'C315', 'C317']],
    'CHL1':[['C3', 'C17']],
    'ch1m':[['C20', 'C12']],
    'CHIM':[['C20', 'C12']],
    'CH0M':[['C20', 'C12']],
    'ERG':[['C3', 'C17']],
    }

CENTRAL_ATOM_OF = {
    'PC':'P',
    'PE':'P',
    'PS':'P',
    'PI':'P',
    'PA':'P',
    'CHL1':'O3',
    'CHIM':'C20',
    'CH0M':'C20',
    'ch1m':'C20',
    'ERG':'O3',
    'WSC1':'N',
    }

# ====================================================
# Protein specific definitions
# ====================================================

RESNAME_SEQUENCE_OF_PROTEIN = {
    "WSC1":[ 'A NVGAI VGGVV GGVVG AVAIA LCILL IVRHI N' ],
    }
