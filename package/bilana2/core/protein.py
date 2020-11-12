'''
    ================================================================================================
    Core objects:  --- :mod: 'bilana.core.protein
        - Protein
            attributes
            functions
        - get_residuegroup_from_seq         Returns group of mda.residues 
                                            for specified 1-lettercoded protein sequence
    
    !!! INCOMPLETE !!!

    ================================================================================================
'''


import os
import logging
import numpy as np
import pandas as pd
import MDAnalysis as mda
import MDAnalysis.analysis.helanal as hl

from bisect import bisect

from .neighbor import get_ref_positions

LOGGER = logging.getLogger("bilana2.core.protein")

AMINO_ACIDS = ["ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLU", "GLN", "GLX", "GLY", "HSD",
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
            "HSD":"H",
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


def get_residuegroup_from_seq(protein_sequences, universe: mda.Universe) -> mda.ResidueGroup:
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
    LOGGER.debug("protein selection is %s", protein_selection)

    ### Translate 3 letter code to 1 letter sequence ###
    all_sequence_1letter = ""
    for resname in protein_selection.resnames:
        all_sequence_1letter += LETTERCODE_3to1[resname]
    LOGGER.debug("protein sequence(s) found: %s", all_sequence_1letter)

    for sequence in sorted(protein_sequences, key=len, reverse=True):
        ### Look for specified <protein_sequence> in whole protein sequence ###
        search = 1
        LOGGER.debug("looking for protein sequence %s", sequence)
        while search:
            ind = all_sequence_1letter.find(sequence)

            if ind >= 0:
                LOGGER.debug("Appending residue group: %s",  protein_selection[ ind:ind+len(sequence) ])
                residue_groups.append( protein_selection[ ind:ind+len(sequence) ] )

                protein_selection = protein_selection[:ind] + protein_selection[ind+len(sequence):]
                all_sequence_1letter = all_sequence_1letter[:ind] + all_sequence_1letter[ind+len(sequence):]
                LOGGER.debug("protein_selection now is %s and lettercode %s", protein_selection, all_sequence_1letter )

                search += 1

            else:
                if search == 1:
                    raise LookupError("input sequence {} not found in bilayer".format(sequence))
                search = False

    if len(protein_selection) != 0:
        oneletterseq = ''.join([LETTERCODE_3to1[resname] for resname in protein_selection.resnames])
        raise ValueError("There were protein residues found which could not be assigned to a sequence:\n{}\nsequence: {}".format(protein_selection, oneletterseq))

    LOGGER.debug("returning residue_groups %s", residue_groups)
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
            universe,
            outputfilepath="./",
            selection="name CA",
            origin_pdbfile="origin.pdb",
            matrix_filename="bending_matrix.dat",
            summary_filename="summary.txt",
            screw_filename="screw.xvg",
            tilt_filename="local_tilt.xvg",
            bend_filename="local_bend.xvg",
            twist_filename="unit_twist.xvg",
            fitted_tilt_filename="fit_tilt.xvg",
            ref_axis=None,
            ):
        origin_pdbfile   = origin_pdbfile.replace(".pdb", str(self.id)+".pdb")
        matrix_filename  = matrix_filename.replace(".dat", str(self.id)+".dat")
        summary_filename = summary_filename.replace(".txt", str(self.id)+".txt")
        screw_filename   = screw_filename.replace(".xvg", str(self.id)+".xvg")
        tilt_filename    = tilt_filename.replace(".xvg", str(self.id)+".xvg")
        bend_filename    = bend_filename.replace(".xvg", str(self.id)+".xvg")
        twist_filename   = twist_filename.replace(".xvg", str(self.id)+".xvg")
        fitted_tilt_filename = fitted_tilt_filename.replace(".xvg", str(self.id)+".xvg")

        path_to_analysis_files = "{}helanal/".format(outputfilepath)
        os.makedirs(path_to_analysis_files, exist_ok=True)
        prefix = path_to_analysis_files

        hl.helanal_trajectory(universe,
            select=selection,
            residuegroup=self.residues,
            origin_pdbfile=origin_pdbfile,
            matrix_filename=matrix_filename,
            summary_filename=summary_filename,
            screw_filename=screw_filename,
            tilt_filename=tilt_filename,
            bend_filename=bend_filename,
            twist_filename=twist_filename,
            fitted_tilt_filename=fitted_tilt_filename,
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
        LOGGER.debug("leafletregions %s", leafletregions)
        for atm in self.reference_atoms:
            pos = atm.positions[0]
            region_ind = bisect(leafletregions, pos[2])
            leaflet = 0 if region_ind < 2 else 1
            resid_to_leaflet[atm.resids[0]] = leaflet
        LOGGER.debug("resid to leaflet %s", resid_to_leaflet)
        return resid_to_leaflet


def create_protein_neighborfile(systeminfo,
    rcut=14,
    outputfilename="protein_neighbors.dat",
    overwrite=False):
    '''
        rcut is in Angstrom!
        Calculate neighborhood of protein:
            protid   time     Number_of_neighbors      List_of_Neighbors
    '''

    refatoms = systeminfo.reference_atomselection 

    lipid_resnames = systeminfo.pl_resnames + systeminfo.sterol_resnames
    PL_refatomnames   = ' '.join([systeminfo.ff.central_atom_of(resname) \
        for resname in lipid_resnames])
    prot_refatomnames = systeminfo.ff.central_atom_of("protein")

    for protein in systeminfo.proteins:
        protid = protein.id
        LOGGER.info("At protein %s", protid)

        outputfilename = outputfilename.replace(".dat", "_{}.dat".format(protid))

        traj_len = len(systeminfo.universe.trajectory)

        if not overwrite and os.path.isfile(outputfilename):
            LOGGER.info("File {} exists, will not overwrite".format(outputfilename))
            return

        with open(outputfilename, "w") as outf:
            print("{: <20}{: <20}{: <20}{: <20}{: <20}".format("resid", "time", "leaflet", "Number_of_neighbors", "List_of_Neighbors"), file=outf)
            for t in range(traj_len):
                time = systeminfo.universe.trajectory[t].time
                if not systeminfo.within_timerange(time):
                    continue
                LOGGER.info("At time %s", time)
                refatomgrp = systeminfo.universe.select_atoms(refatoms)
                refpositions = get_ref_positions(systeminfo, "atom", refatomgrp) # leaflets=[(resid1, pos1), ...]

                LOGGER.debug("Found %s atoms: %s ", len(refatomgrp.atoms), refatomgrp.atoms)
                LOGGER.debug("Dimension of leaflets %s", np.array(refpositions).shape)

                position_array   =  np.array([pos for resid, pos in refpositions])
                position_array2d =  position_array.copy()
                position_array2d[:,2] = 0



                for refatm in protein.reference_atoms:

                    host_pos = refatm.positions[0]
                    host_pos2d = host_pos.copy()
                    host_pos2d[2] = 0
                    hostid = refatm.resids[0]

                    dist_array = mda.lib.distances.distance_array(
                        host_pos,
                        position_array,
                        box=systeminfo.universe.dimensions)[0]
                    dist_array2d = mda.lib.distances.distance_array(
                        host_pos2d,
                        position_array2d,
                        box=systeminfo.universe.dimensions)[0]

                    leaflet = protein.resid_to_leaflet[hostid]

                    #neibresidues = systeminfo.universe.atoms.select_atoms(
                    #        "name {} and around {} (name {} and resid {})"\
                    #        .format(PL_refatomnames, rcut, prot_refatomnames, resid)).residues

                    neiblist = []
                    for resndx, distance in enumerate(dist_array):
                        neibid = refpositions[resndx][0]
                        neib_leaf = systeminfo.convert.resid_to_leaflet[neibid]

                        if leaflet == neib_leaf:
                            distance2d = dist_array2d[resndx]
                            if distance2d <= rcut:
                                neiblist.append( neibid )

                    neiblist = list(set(neiblist)) # delete duplicates
                    neiblist.sort()
                    neiblist = [ str(n) for n in neiblist if n != hostid ] # delete host entry
                    n_neibs = len(neiblist)
                    neiblist = ','.join( neiblist )
                    LOGGER.debug("Neiblist: %s", neiblist)


                    print("{: <20}{: <20}{: <20}{: <20}{: <20}".format(hostid, time, leaflet, n_neibs, neiblist), file=outf)



def create_protein_protein_distancefile(systeminfo, protein1, protein2, outputfilename="tmd_distance.csv"):
    ''' '''
    u = systeminfo.universe

    outputfilename = outputfilename.replace(".csv", "{}-{}.csv".format(protein1.id, protein2.id))
    
    times = []
    distances = []
    for t in range(len(u.trajectory)):
        time = systeminfo.universe.trajectory[t].time
        if not systeminfo.within_timerange(time):
            continue
        LOGGER.info("at time %s", time)

        # take COM as reference positions
        com1 = protein1.residues.atoms.center_of_mass()
        com2 = protein2.residues.atoms.center_of_mass()

        # calculate distance
        dist = np.linalg.norm(com1 - com2) / 10 ## ouput in nm

        # output data
        times.append(time)
        distances.append(dist)
    dat = pd.DataFrame({"time":times, "distance":distances})
    dat.to_csv(outputfilename, index=False)




