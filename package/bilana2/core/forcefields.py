"""
    ===========================================================
    Core object: Forcefield  --- :mod: 'bilana.core.forcefields
    ===========================================================

    This module stores all information regarding force field atom names defined in respective
    files in ..definitions.
    Provided are functions that easily allow for gathering specific atom groups.


"""

import logging
from ..definitions import charmm
from ..definitions import martini

LOGGER = logging.getLogger("bilana2.core.forcefields")

PROTEIN_SINGLE_LETTER_TO_RESNAME = {}

class Forcefield(object):
    """
        Manager for the different force fields

        Functions to retrieve lists of atoms:
            - central_atom_of( resname )
            - head_atoms_of( resname )
            - tail_atoms_of( resname )
            - tailcarbons_of( resname )
            - tailhydrogens_of( resname )
            - scd_tail_atoms_of( resname )

    """
    def __init__(self, ffname):
        if ffname == "charmm":
            self._ff = charmm
        elif ffname == "martini":
            self._ff = martini
        else:
            raise ValueError("{} does not exist".format(ffname))

    def is_sterol(self, resname):
        return resname in self._ff.STEROLS
    def is_protein(self, resname):
        return resname in self._ff.PROTEIN_SEQUENCES
    def is_water(self, resname):
        return resname in self._ff.WATER

    def central_atom_of(self, resname):
        ''' Return the central atom of specified lipid as string '''
        if self.is_sterol(resname):
            return self._ff.CENTRAL_ATOM_OF[resname]
        elif self.is_protein(resname):
            return 'N'
        else:
            head = resname[-2:]
            return self._ff.CENTRAL_ATOM_OF[head]

    def head_atoms_of(self, resname):
        ''' Returns a list of all head atoms for specified lipid species '''
        if self.is_sterol(resname):
            return self._ff.HEAD_ATOMS_OF[resname]
        else:
            head = resname[-2:]
            return self._ff.HEAD_ATOMS_OF[head]

    def tail_atoms_of(self, resname):
        ''' Returns a list of all tail atoms for specified lipid species '''
        if self.is_sterol(resname):
            LOGGER.warning("WARNING: No hydrogens are included")
            return self._ff.TAIL_ATOMS_OF[resname]
        else:
            tail = resname[:-2]
            return self._ff.TAIL_ATOMS_OF[tail]

    def tailcarbons_of(self, resname):
        ''' Returns only carbon atoms of tail '''
        if self.is_sterol(resname):
            return self._ff.TAIL_ATOMS_OF[resname]
        else:
            tail = resname[:-2]
            return self._ff.TAILCARBONS[tail]

    def tailhydrogens_of(self, resname):
        ''' Returns only hydrogen atoms of tail '''
        if self.is_sterol(resname):
            raise ValueError("Hydrogens of sterol not (yet?) included.")
        else:
            tail = resname[:-2]
            return self._ff.TAILHYDROGENS[tail]

    def scd_tail_atoms_of(self, resname):
        ''' Returns a list of relevant carbons for calculation of scd '''
        if self.is_sterol(resname):
            return self._ff.SCD_TAIL_ATOMS_OF[resname]
        elif self.is_protein(resname):
            return [['N', 'C']]
        else:
            tail = resname[:-2]
            return self._ff.SCD_TAIL_ATOMS_OF[tail]

    def get_resnames_from_sequence(self, sequence):
        ''' returns a list of residue names from sequence '''
        raise NotImplementedError("Function not yet implemented")

    def _check_forcefield_file(self):
        ''' Check wether all components were added properly'''
