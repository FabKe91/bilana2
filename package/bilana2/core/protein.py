
'''

'''
import os
import logging
import MDAnalysis as mda
import MDAnalysis.analysis.helanal as hl

from bisect import bisect

LOGGER = logging.getLogger("bilana2.core.protein")

AMINO_ACIDS = ["ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ]

LETTERCODE_3to1 = {
            "ALA":"A",
            "ARG":"R",
            "ASN":"N",
            "ASP":"D",
            "ASX":"B",
            "CYS":"C",
            "GLU":"E",
            "GLN":"Q",
            "GLX":"Z",
            "GLY":"G",
            "HIS":"H",
            "ILE":"I",
            "LEU":"L",
            "LYS":"K",
            "MET":"M",
            "PHE":"F",
            "PRO":"P",
            "SER":"S",
            "THR":"T",
            "TRP":"W",
            "TYR":"Y",
            "VAL":"V",
     }
LETTERCODE_1to3 = dict( reversed(item) for item in LETTERCODE_3to1.items() )


def get_residuegroup_from_seq(protein_sequences, universe: mda.Universe) -> mda.residues:
    ''' Returns group of mda.residues for specified 1-lettercoded protein sequence
        input:
            protein_sequences: 1-letter-coded AA sequence that are found in system
        ToDo:
            If two sequences with fitting beginning are in system:
                - start with longer sequence
                - look for shorter sequence with longer one deleted

    '''
    residue_groups = []

    protein_selection = universe.atoms.select_atoms("protein").residues

    ### Translate 3 letter code to 1 letter sequence ###
    all_sequence_1letter = ""
    for resname in protein_selection.resnames:
        all_sequence_1letter += LETTERCODE_3to1[resname]

    for sequence in sorted(protein_sequences, key=len, reverse=True):
        ### Look for specified <protein_sequence> in whole protein sequence ###
        search = 1
        while search:
            ind = all_sequence_1letter.find(sequence)

            if ind >= 0:
                residue_groups.append( protein_selection[ ind:ind+len(sequence) ] )

                protein_selection = protein_selection[:ind] + protein_selection[ind+len(sequence)]
                all_sequence_1letter = all_sequence_1letter[:ind] + all_sequence_1letter[:ind+len(sequence)]

                search += 1

            else:
                if search == 1:
                    raise LookupError("input sequence {} not found in bilayer".format(sequence))
                search = False

    if len(protein_selection) != 0:
        raise ValueError("There were protein residues found which could not be assigned to a sequence")

    return residue_groups



class Protein():
    '''
        Class Protein holds together a group of mda.residues of a specified sequence
        Input:


        Attributes:

    '''
    def __init__(self, prot_id, residues, forcefield, bilayerregions):
        self.id       = prot_id
        self.residues = residues
        self.sequence = self._sequence_from_residues()
        self.resids   = residues.resids
        self.resnames = residues.resnames
        self.ff       = forcefield

        self.reference_atoms  = self._add_reference_atoms()
        self.resid_to_leaflet = self._assign_resids_to_leaflet(bilayerregions)

    def run_helanal(self,
            outputfilepath="./",
            selection="name CA",
            origin_pdbfile="origin.pdb",
            matrix_filename="bending_matrix.dat",
            summary_filename="summary.txt",
            screw_filename="screw.xvg",
            tilt_filename="local_tilt.xvg",
            bend_filename="local_bend.xvg",
            twist_filename="unit_twist.xvg",
            ref_axis=None,
            ):
        path_to_analysis_files = "{}tmdanalysis/".format(outputfilepath)
        os.makedirs(path_to_analysis_files, exist_ok=True)
        prefix = path_to_analysis_files
        hl.helanal_trajectory(self.residues, ### !!! Dont know if that works !!!
            origin_pdbfile=origin_pdbfile,
            matrix_filename=matrix_filename,
            summary_filename=summary_filename,
            screw_filename=screw_filename,
            tilt_filename=tilt_filename,
            bend_filename=bend_filename,
            twist_filename=twist_filename,
            ref_axis=ref_axis,
            prefix=prefix,
            )

    #def _define_helanal_outputfilesnames(self): ### !!! UNDER CONSTRUCTION !!!
    #    ''' '''
    #    self.selection        = "name CA"
    #    self.origin_pdbfile   = "origin.pdb"
    #    self.matrix_filename  = "bending_matrix.dat"
    #    self.summary_filename = "summary.txt"
    #    self.screw_filename   = "screw.xvg"
    #    self.tilt_filename    = "local_tilt.xvg"
    #    self.bend_filename    = "local_bend.xvg"
    #    self.twist_filename   = "unit_twist.xvg"


    def _sequence_from_residues(self):
        ''' returns protein 1-letter-coded sequence derived from residuegroup '''
        sequence = ""
        for resname in self.residues.resnames:
            sequence += LETTERCODE_3to1[resname]
        return sequence

    def _add_reference_atoms(self):
        ''' Collects all reference atoms defined in ff and stores them in a list '''
        refatms = []
        for residue in self.residues:
            atmsel = residue.atoms.select_atoms("name {}".format(self.ff.central_atom_of(residue.resname)))
            refatms.append(atmsel)
        return refatms
    def _assign_resids_to_leaflet(self, leafletregions):
        ''' Use definition of leafletregions=[ h1/t1, t1/t2, t2/h2, ] to assign prot resids '''
        resid_to_leaflet = {}
        for atm in self.reference_atoms:
            pos = atm.positions
            region_ind = bisect(leafletregions, pos[2])
            leaflet = 0 if region_ind < 2 else 1
            resid_to_leaflet[atm.resid] = leaflet
        return resid_to_leaflet


def create_protein_neighborfile(systeminfo,
    rcut=1.4,
    outputfilename="protein_neighbors.csv",
    overwrite=False):
    '''
        Calculate neighborhood of protein:
            protid   time     Number_of_neighbors      List_of_Neighbors
    '''

    lipid_resnames = systeminfo.PL_resnames + systeminfo.sterol_resnames
    PL_refatomnames   = ' '.join([systeminfo.ff.central_atom_of(resname) \
        for resname in lipid_resnames])
    prot_refatomnames = systeminfo.ff.central_atom_of("protein")

    for protein in systeminfo.proteins:
        protid = protein.id

        traj_len = len(systeminfo.universe.trajectory)

        if not overwrite and os.path.isfile(outputfilename):
            LOGGER.info("File {} exists, will not overwrite".format(outputfilename))
            return

        with open(outputfilename, "w") as outf:
            print("{: <20}{: <20}{: <20}{: <20}".format("resid", "time", "Number_of_neighbors", "List_of_Neighbors"), file=outf)
            for t in range(traj_len):
                time = systeminfo.universe.trajectory[t].time
                if not systeminfo.within_timerange(time):
                    continue
                neibresidues = systeminfo.universe.atoms.select_atoms(
                        "name {} and around {} name {}"\
                        .format(PL_refatomnames, rcut, prot_refatomnames)).residues

                neiblist = ','.join( [ str(resid) for resid in neibresidues.resids ] )
                n_neibs  = len(neibresidues)
                print("{: <20}{: <20}{: <20}{: <20}".format(protid, time, n_neibs, neiblist), file=outf)



def create_protein_protein_distancefile(protein1, protein2):
    ''' '''
