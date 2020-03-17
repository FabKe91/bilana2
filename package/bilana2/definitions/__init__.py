'''
    definitions gathers all force field information

    Obligatory definitions per force field are:
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
        - PROTEIN_RESIDUES      Lists all known protein residue types #### In future all residue
                                types should be included
        - SOLVENTS              Lists the names of all solvent types in the respective force field

'''
from . import charmm
from . import martini
