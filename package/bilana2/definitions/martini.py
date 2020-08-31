
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

TAILS             = ['DP', 'DI', 'DU']
HEADS             = ['PC']
STEROLS           = ['CHL1', 'ERG']
PROTEIN_SEQUENCES = []
AMINO_ACIDS       = ['VAL', 'GLY', 'ALA', 'ILE', 'LEU', 'CYS', 'ARG', 'HSD']
IONS              = ["CL", "POT", "NA" ]
WATER             = ["TIP3", "SOL"]
SOLVENTS          = WATER + IONS

# =================================================
# Glycerolpart, also carbonylpart of FA is included
# =================================================

GLYCEROLATOMS = ['GL1', 'GL2', ]

# =================================================
# Head atom definition
# =================================================

HEAD_ATOMS_OF = {
    'PC':['PO4', 'NC3'],
    'CHOL':['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'],
    'ERG':['ROH', 'R1', 'R2', 'R3', 'R4', 'R5'],
    }

# =================================================
# Tail definitions
# =================================================

TAILCARBONS = {
    'DP':[['C1A', 'C2A', 'C3A', 'C4A'],                 #16:0  ### This length should be a range as martini cannot distinguish between particular lengths
          ['C1B', 'C2B', 'C3B', 'C4B']],                #16:0  ### This length should be a range as martini cannot distinguish between particular lengths
    'DI':[['C1A', 'D2A', 'D3A', 'C4A'],                 #14:0  ### This length should be a range as martini cannot distinguish between particular lengths
          ['C1B', 'D2B', 'D3B', 'C4B']],                #14:0  ### This length should be a range as martini cannot distinguish between particular lengths
    'DU':[['C1A', 'D2A', 'D3A', 'C4A'],                 #14:0  ### This length should be a range as martini cannot distinguish between particular lengths
          ['C1B', 'D2B', 'D3B', 'C4B']],                #14:0  ### This length should be a range as martini cannot distinguish between particular lengths
    'CHL1':[['C1', 'C2']],
    'ERG':[['C1', 'C2']],
    }

TAILHYDROGENS = [] # Martini only includes heavy atoms

TAIL_ATOMS_OF = TAILCARBONS

# ====================================================
# In the following special groups of atoms are defined
# ====================================================

SCD_TAIL_ATOMS_OF = {
    'DP':[TAILCARBONS['DP'][0][::2], TAILCARBONS['DP'][1][::2]],
    'DI':[TAILCARBONS['DI'][0][::2], TAILCARBONS['DI'][1][::2]],
    'DU':[TAILCARBONS['DI'][0][::2], TAILCARBONS['DI'][1][::2]],
    'CHL1':[['ROH', 'C2']],
    'ERG':[['ROH', 'C2']],
    }

CENTRAL_ATOM_OF = {
    'PC':'PO4',
    'CHL1':'ROH',
    'ERG':'ROH',
    }

# ====================================================
# Protein specific definitions
# ====================================================

RESNAME_SEQUENCE_OF_PROTEIN = {}
