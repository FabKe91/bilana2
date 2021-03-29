'''
==========================================================
Core object: Systeminfo  --- :mod: 'bilana.core.systeminfo
==========================================================

'''
import os
import re
import logging

import numpy as np
import MDAnalysis as mda

from .protein import Protein, get_residuegroup_from_seq
from ..lib import common as cm
from .forcefields import Forcefield


LOGGER = logging.getLogger("bilana2.core.systeminfo")


class Systeminfo(object):
    """
    The bilana Systeminfo class stores all information of a bilayer simulation using
    Gromacs MD package
    It is supposed in aiding common tasks like
        - easy handling of all mdfiles for gromacs
          analysis application (all paths stored in Systeminfo.files attribute)
        - gathering all force field information like atom names and atoms used for
          order parameter calculation

    Systeminfo class instances are supposed to be the main input for all following analysis like
    - order parameter calculation (see order.py)
    - mapping the neighborhood of all lipids (see neighbor.py)
    - lipid interaction energy (see energy.py)

    Currently only two force-fields are included:
        - CHARMM36
        - Martini

    Example layout for input file:

        #system info
        system: dppc_chol20
        lipidmolecules: DPPC,CHL1
        temperature: 290
        timeframe: 0,900000000,1000

        cutoff: 1.0
        forcefield: charmm36

        refatomselection: (resname DPPC DUPC and name P) or (resname CHL1 and name O3) or (resname CHIM and name C20)

        # absolute path to mdfiles
        mdfiles: /scratch/f_kell07/mdfiles



    """

    def __init__(self, inputfilepath="inputfile"):
        # ============================================================
        # Following will be set in self.read_input_and_set_attributes
        # ============================================================

        self.system                  = None
        self.temperature             = None
        self.cutoff                  = None
        self.reference_atomselection = None
        self.molecules               = None
        self.forcefieldname          = None
        self.mdfilespath             = None
        self._inputfile_times        = None
        self.protein_sequences       = None

        self._read_input_and_set_attributes( inputfilepath )

        # ==============================================================
        # Setting more complex attributes or initialize class instances.
        # For further information look into docs of respective classes
        # ==============================================================

        self.path     = Files(inputfilepath, self.mdfilespath, self.system, self.temperature)
        self.ff       = Forcefield(self.forcefieldname)
        self.universe = mda.Universe(self.path.gro, self.path.trj)
        self.convert  = Conv(self.universe, self.ff, self.molecules)

        # ===========================================
        # Store geometric information
        # ===========================================

        self.bilayerregions = self._define_bilayer_regions()

        # ==============================================================
        #  Set the times from input file or gathered from mda.Universe
        # ==============================================================

        self.t_start  = None
        self.t_end    = None
        self.dt       = None

        self._set_times(self.universe)

        # =====================================
        # Store system composition information
        # =====================================

        self.lipid_resids     = np.unique( self.universe.select_atoms("resname {}".format( ' '.join( self.molecules ) ) ).resids )
        self.number_of_lipids = self.lipid_resids.size
        self.solvent_resnames = [mol for mol in np.unique(self.universe.atoms.resnames) if mol not in self.molecules]
        self.water_resname    = [mol for mol in np.unique(self.universe.atoms.resnames) if self.ff.is_water(mol)][0]

        self.pl_resnames      = [ mol for mol in np.unique(self.universe.residues.resnames) if \
                    ( (mol in self.molecules) and not ( self.ff.is_protein(mol) or self.ff.is_sterol(mol) ) ) ]
        self.pl_resids        = np.unique(self.universe.atoms.select_atoms( 'resname {}'.format( ' '.join( self.pl_resnames ) ) ).resids)

        self.sterol_resnames  = [ mol for mol in np.unique(self.universe.residues.resnames) if \
                    ( (mol in self.molecules) and self.ff.is_sterol(mol) ) ]
        self.sterol_resids    = []

        if self.sterol_resnames:
            self._add_sterol_resids()

        self._check_if_all_molecules_found()

        # ===========================================
        # Store information about proteins in system
        # ===========================================

        self.proteins = []
        self.n_proteins = len(self.protein_sequences)

        if self.protein_sequences:
            self._add_proteins()


    def within_timerange(self, time):
        return (self.t_start <= time <= self.t_end) and (time % self.dt == 0)

    def _check_if_all_molecules_found(self):
        found_lipids = self.pl_resnames + self.sterol_resnames + self.protein_resnames
        not_found    = list(set(self.molecules) - set(found_lipids))
        if not_found:
            raise ValueError("Not all lipids given in inputfile found in structure file! Incorrect name?\n"\
                    "Lipid(s) not found: {}".format(' '.join(not_found) ) )

    def _add_sterol_resids(self):
        self.sterol_resids += list(np.unique(self.universe.atoms.select_atoms(
            'resname {}'.format( ' '.join( self.sterol_resnames ) )
            ).resids))

    def _add_proteins(self):
        for prot_id, seq in enumerate(self.protein_sequences):
            protein_residues = get_residuegroup_from_seq(seq, self.universe)
            self.proteins.append( Protein(prot_id, protein_residues, self.ff, self.bilayerregions) )

    def _read_input_and_set_attributes(self, inputfilepath ):
        ''' Reading input file
                - parameters are separated by a colon
                - comments can be made using #
                - multiple input parameters are separated by a comma
                - refatomselection uses MDAnalysis selection syntax
        '''
        systeminfo = {}
        with open(inputfilepath,"r") as inputf:
            # Creates a list like [[system,dppc_chol],[temperature,290]]
            regex = re.compile(r'^([\w,\d,\s]*):([\w, \d, \s, \(, \), \., /]*)\s*#*.*$')
            for line in inputf:
                match = regex.match(line)
                if match is not None:
                    key = match.group(1).strip()
                    item = match.group(2).strip().replace("\n", "")
                    if key != "refatomselection":
                        item = item.replace(" ", "")
                    systeminfo[key] = item

        self.system                  = systeminfo.pop("system")
        self.temperature             = int(systeminfo.pop("temperature"))
        self.cutoff                  = float(systeminfo.pop("cutoff"))
        self.reference_atomselection = systeminfo.pop("refatomselection")
        self.molecules               = [i.strip() for i in systeminfo.pop("lipidmolecules").split(",")]
        self.mdfilespath             = systeminfo.pop("mdfiles")
        self.forcefieldname          = systeminfo.pop("forcefield")
        self._inputfile_times        = systeminfo.pop("timeframe").split(",")
        self.protein_sequences       = systeminfo.pop("protein_sequences").split(",")

        self._check_inputfile(systeminfo)

    def _check_inputfile(self, remaining_input):
        ''' Check wether obligatory input parameters are missing
            or wrong additional parameters were parsed
        '''
        obligatory_input_parameters = {
            "system":self.system,
            "temperature":self.temperature,
            "cutoff":self.cutoff,
            "refatomselection":self.reference_atomselection,
            "lipidmolecules":self.molecules,
            "forcefield":self.forcefieldname,
            "mdfiles":self.mdfilespath,
            "timeframe":self._inputfile_times,
        }

        # Check wether key attributes remain unset
        if None in obligatory_input_parameters.values():
            missing_parameters = []
            for key, val in obligatory_input_parameters.items():
                if val is None:
                    missing_parameters.append(key)
            raise AttributeError( "Following parameters are missing:\n{}"\
                .format( '\n'.join(missing_parameters) ) )

        # Raises an exception if additional/unknown parameters are parsed
        if remaining_input:
            raise AttributeError( "Following parameters remain unused. Typo?\n{}"\
                .format( '\n'.join(remaining_input) ) )

    def _set_times(self, universe):
        '''
        Sets all time attributes:
            self.dt
            self.t_end
            self.t_start
        '''
        t_end_mda   = int( universe.trajectory[-1].time )
        dt_mda      = int( universe.trajectory.dt )
        t_end_inp   = int( self._inputfile_times[1] )
        dt_inp      = int( self._inputfile_times[2] )

        self.t_start = int( self._inputfile_times[0] )
        self.dt      = max(dt_inp, dt_mda)
        self.t_end   = min(t_end_inp, t_end_mda)

    def _define_bilayer_regions(self):
        ''' Define borders of regions within bilayer
            like [ h1/t1, t1/t2, t2/h2, ] from first frame
            assumed lipid proportions are

              P       |   Head region: 1/4
              |_      |_____ Total: 2.5 nm
             /  \     |  Tail region: 3/4
            /    \    |

        '''
        pl_resnames = ' '.join(self.pl_resnames)
        pl_refatomnames = ' '.join(
            [self.ff.central_atom_of[resname] for resname in self.pl_resnames])
        pl_refatoms = self.universe.atoms.select_atoms(
            "resname {} and name {}".format(pl_resnames, pl_refatomnames))
        zpos_refatms = pl_refatoms.positions[:,2]

        bilayercenter = zpos_refatms.mean()
        upperavg = zpos_refatms[zpos_refatms >= bilayercenter]
        loweravg = zpos_refatms[zpos_refatms <  bilayercenter]

        upperborder = (bilayercenter + upperavg * 0.875)  # Tail region + 0.5*headregion
        lowerborder = (bilayercenter - loweravg * 0.875)  #  3/4        +   0.5 * 1/4

        return [lowerborder, bilayercenter, upperborder]



class Conv(object):
    """
        Stores all conversion dicts
    """
    def __init__(self, universe: mda.Universe, ff: Forcefield, molecules: list):
        self.index_to_resid    = np.array([0, *universe.atoms.resids])
        self.resid_to_resname  = self._get_resid_to_resname(universe)
        self.resid_to_leaflet  = self._leaflet_assignment_of_frame(
                                    universe, ff, molecules)
        self.leaflet_to_resids = [[], []]
        for key, val in self.resid_to_leaflet.items():
            self.leaflet_to_resids[val].append(key)

    def _get_resid_to_resname(self, u: mda.Universe) -> dict:
        ''' Create a dictionary that maps all resids to the respective molecules resname '''
        mask     = [True, *np.array(np.diff(u.atoms.resids), dtype=bool)]
        resids   = u.atoms.resids[ mask ]
        resnames = u.atoms.resnames[ mask ]
        return dict( zip( resids, resnames ) )

    def _leaflet_assignment_of_frame(self,
        universe: mda.Universe, ff, molecules) -> dict:
        '''
            Input: MDAnalysis.AtomGroup with selected atoms to assign a leaflet to
            Returns a pandas dataframe with entries
                <resid> <resname> <leaflet>
        '''
        resids = []
        leaflet_assignment = []

        for residue in universe.residues:
            if residue.resname not in molecules:
                continue

            headnames = ' '.join( ff.head_atoms_of(residue.resname) )
            ### unpacking ###
            tailcarbonlist = [i for j in ff.tailcarbons_of(residue.resname) for i in j]
            tailnames = ' '.join( tailcarbonlist )

            coord_head = residue.atoms.select_atoms( "name {}".format( headnames )
                    ).atoms.center_of_geometry()
            coord_tail = residue.atoms.select_atoms( "name {}".format( tailnames )
                    ).atoms[:10].center_of_geometry()

            leaflet = cm.molecule_leaflet_orientation(coord_head, coord_tail)
            resids.append(residue.resid)
            leaflet_assignment.append(leaflet)

        return dict( zip(resids, leaflet_assignment ) )

class Files(object):
    """
        Stores paths to MD files form Gromacs used in further analysis

        Files must be stored in the mdfiles folder given in the inputfile in the following manner:
            mdfiles/
                |--md_trj/*_whole.xtc
                |--initial_coords/*.gro
                |--psf/*.top
                |--tpr/*.tpr
                |--enr/*.edr

        MD trajectories must be made whole and in .xtc format.


    """

    def __init__(self, inputfilepath, mdfilespath, system, temperature):
        inpdir = os.path.dirname(inputfilepath)
        if not inpdir:
            inpdir = '.'

        self.mdfiles = mdfilespath

        self.trj    = '{}/md_trj/{}_{}_whole.xtc'.format(mdfilespath, system, temperature)
        self.gro    = '{}/initial_coords/{}.gro'.format(mdfilespath, system)
        self.top    = '{}/psf/{}.top'.format(mdfilespath, system)
        self.tpr    = '{}/tpr/{}_{}.tpr'.format(mdfilespath, system, temperature)
        self.edr    = '{}/enr/{}_{}.edr'.format(mdfilespath, system, temperature)

        self.index  = "{}/indexfiles".format(inpdir)
        self.data   = "{}/datafiles".format(inpdir)
        self.tmp    = "{}/tempfiles".format(inpdir)
        self.energy = "{}/energyfiles".format(inpdir)
        self.log    = "{}/logfiles".format(inpdir)

        self._check_exist_files()

    def create_folders(self):
        '''
            Creates following folders:
                - datafiles     datafiles stores all analysis data
                - indexfiles    indexfiles stores gromacs .ndx files used for the energy calculation
                - tempfiles     tempfiles stores all temporary files
                - energyfiles   energyfiles stores all files used in energy calculations:
                                *.mdp, *.tpr, *.edr, *.xvg
        '''
        os.makedirs(self.index, exist_ok=True)
        os.makedirs(self.data, exist_ok=True)
        os.makedirs(self.tmp, exist_ok=True)
        os.makedirs(self.energy, exist_ok=True)
        os.makedirs(self.log, exist_ok=True)

    def _check_exist_files(self):
        def exists(fname):
            return os.path.isfile(fname)

        fnames = [self.gro, self.tpr]
        missing_files = []
        trj_exists    = exists( self.trj )

        for fname in fnames:
            if not exists(fname):
                missing_files.append(fname)

        err_str = ""
        if missing_files:
            err_str += "Following files are missing: {}\n".format( ' '.join(missing_files) )
        if not trj_exists:
            err_str += "Trajectory {} is missing. Are you sure you provided"\
                " a trajectory made whole?".format(self.trj)
        if err_str:
            raise FileNotFoundError(err_str)
