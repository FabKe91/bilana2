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

from MDAnalysis.lib.distances import distance_array, calc_bonds, self_distance_array, calc_angles

from bisect import bisect

from .neighbor import get_ref_positions

LOGGER = logging.getLogger("bilana2.core.protein")
#LOGGER.setLevel("DEBUG")

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
    LOGGER.debug("protein selection is %s with %s items", protein_selection, len(protein_selection))

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
                protein_group = protein_selection[ ind:ind+len(sequence) ]
                LOGGER.debug("Appending residue group: %s",  protein_group)
                LOGGER.debug("Atom ids %s", protein_group.ids)
                residue_groups.append( protein_group )
                LOGGER.debug("Having now %s residue groups", len(residue_groups))

                protein_selection = protein_selection[:ind] + protein_selection[ind+len(sequence):]
                all_sequence_1letter = all_sequence_1letter[:ind] + all_sequence_1letter[ind+len(sequence):]
                LOGGER.debug("protein_selection now is %s (%s items) and lettercode %s", protein_selection, len(protein_selection),  all_sequence_1letter )

                search += 1

                LOGGER.debug("At search %s", search)

            else:
                if search == 1:
                    print(search, "before error")
                    raise LookupError("input sequence {} not found in bilayer".format(sequence))
                search = False

            if search > 1:
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
            fitted_tilt_vector_filename="fit_tilt_vector.xvg",
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
        fitted_tilt_vector_filename    = fitted_tilt_vector_filename.replace(".xvg", str(self.id)+".xvg")

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
            fitted_tilt_vector_filename=fitted_tilt_vector_filename,
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
    cutoff=4,
    outputfilename="protein_neighbors.dat",
    overwrite=False):
    '''
        cutoff is in Angstrom!
        Calculate neighborhood of protein:
            protid   time     Number_of_neighbors      List_of_Neighbors
    '''

    #refatoms = systeminfo.reference_atomselection 
    #lipid_resnames = systeminfo.pl_resnames + systeminfo.sterol_resnames
    #PL_refatomnames   = ' '.join([systeminfo.ff.central_atom_of(resname) \
    #    for resname in lipid_resnames])
    #prot_refatomnames = systeminfo.ff.central_atom_of("protein")

    lipids = systeminfo.universe.atoms.select_atoms("resname {}".format(' '.join(systeminfo.molecules)))

    if not overwrite and os.path.isfile(outputfilename):
        LOGGER.info("File {} exists, will not overwrite".format(outputfilename))
        return

    with open(outputfilename, "w") as outf:
        print("{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}".format("protid", "resid", "time", "leaflet", "Number_of_neighbors", "List_of_Neighbors"), file=outf)
        for prot in systeminfo.proteins:
            protid = prot.id
            LOGGER.info("At protein %s", protid)

            #refatomgrp = systeminfo.universe.select_atoms(refatoms)
            #refpositions = get_ref_positions(systeminfo, "atom", refatomgrp) # leaflets=[(resid1, pos1), ...]

            for ts in systeminfo.universe.trajectory:
                time = ts.time
                box = ts.dimensions

                if not systeminfo.within_timerange(time):
                    continue
                LOGGER.info("At time %s", time)

                #LOGGER.debug("Found %s atoms: %s ", len(refatomgrp.atoms), refatomgrp.atoms)
                #LOGGER.debug("Dimension of leaflets %s", np.array(refpositions).shape)

                dist_ar = distance_array(prot.residues.atoms.positions, lipids.atoms.positions, box=box)

                ### Get all lipids within CUTOFF for each res
                ndxa = 0
                for i, res in enumerate(prot.residues):

                    leaflet  = prot.resid_to_leaflet[res.resid]

                    ndxb = ndxa + len(res.atoms)

                    #print("from", ndxa, ndxb, dist_ar.shape)

                    atoms_in_range = []
                    for ar in dist_ar[ndxa:ndxb]: # Loop over each atom in residue <res>
                        atoms_in_range.append(lipids.atoms[ar <= cutoff])
                    lipids_in_range = np.sum(atoms_in_range) # returns 0.0 if atom groups are empty
                    if isinstance(lipids_in_range, mda.core.groups.Atom):
                        lipids_in_range = mda.ResidueGroup([lipids_in_range.residue]) # Returns ResidueGroup with 1 residue
                    elif not isinstance(lipids_in_range, float):
                        lipids_in_range = lipids_in_range.residues # This should be the normal case
                    else:
                        lipids_in_range = atoms_in_range[0].residues # Returns ResidueGroup with 0 residues

                    ndxa += len(res.atoms)

                    ### Count lipid types per res
                    #resn_count = []
                    #for mol in systeminfo.molecules:
                    #    resn_count.append(list(lipids_in_range.resnames).count(mol))
                    #resn_count = np.array(resn_count, dtype=int)

                    ### Create a neiblist 
                    neiblist = [str(resid) for resid in lipids_in_range.resids]
                    #        if systeminfo.convert.resid_to_leaflet[resid] == leaflet]

                    neiblist.sort()
                    n_neibs = len(neiblist)
                    neiblist = ','.join( neiblist )
                    #LOGGER.debug("Neiblist: %s", neiblist)
                    LOGGER.debug("id %s | resid %s | time %s | leaf %s | Ntot %s | list %s", protid, res.resid, time, leaflet, n_neibs, neiblist)

                    print("{: <20}{: <20}{: <20}{: <20}{: <20}{: <20}".format(protid, res.resid, time, leaflet, n_neibs, neiblist), file=outf)



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


def calculate_crossing_angle(systeminfo, protein1, protein2,
        outputfilename='crossing_angle.dat',
        binwidth=10000,
        sel="name CA"):


    u = systeminfo.universe
    start = systeminfo.t_start
    end = systeminfo.t_end
    dt = systeminfo.dt

    ### Define selection for begin and end parts of the TMDs ###
    ### The COG is used as reference pos for selection       ###
    sel_b_k = protein1.residues[0:4].atoms.select_atoms(sel)
    sel_b_j = protein2.residues[0:4].atoms.select_atoms(sel)

    sel_e_k = protein1.residues[-4:].atoms.select_atoms(sel)
    sel_e_j = protein2.residues[-4:].atoms.select_atoms(sel)
        
    dt = u.trajectory.dt
    if binwidth < dt:
        binwidth = dt

    start_frame = int(start / dt)
    end_frame = int(end / dt)
    time_interval = int(binwidth / dt)

    data = pd.DataFrame([])

    for i_ts,ts in enumerate(u.trajectory[start_frame:end_frame:time_interval]):

        time = ts.time
        box = u.dimensions

        box_dim = box[:3]        

        b_k, e_k = sel_b_k.center_of_geometry(), sel_e_k.center_of_geometry()
        b_j, e_j = sel_b_j.center_of_geometry(), sel_e_j.center_of_geometry()

        ### t_j-t_k defines the contact vector ###
        t_k, t_j = minimal_distance(b_k,e_k,b_j,e_j)

        crossing_angle = mda.lib.distances.calc_dihedrals(b_k, t_k, t_j, b_j, box=box) * (180/np.pi)

        ### Shift angles accordingly ###
        if crossing_angle > 90:
            crossing_angle = crossing_angle - 90
        if crossing_angle < -90:
            crossing_angle = crossing_angle + 90
        
        #### Why this correction? ####
        dist_com_uncorrected = distance_array(protein1.residues.atoms.center_of_geometry(), protein2.residues.atoms.center_of_geometry())[0][0]

        diffvector1 = np.subtract(b_k,e_k)
        diffvector2 = np.subtract(b_j,e_j)

        tmd_angle = np.arccos(np.dot(diffvector1, diffvector2)/(np.linalg.norm(diffvector1)*np.linalg.norm(diffvector2))) * (180/np.pi)

        if (dist_com_uncorrected > box_dim[0]/2) or (dist_com_uncorrected > box_dim[1]/2):

            s = np.sign(crossing_angle)

            if s != 0:
                crossing_angle = np.abs(tmd_angle)*s

        ### Append data for this timeframe ###
        data = data.append(pd.DataFrame({'time': time, 'crossing_angle': crossing_angle}, index=[i_ts]), ignore_index=False)        
            
    data.to_csv(outputfilename, index=False)


def minimal_distance(begA, endA, begB, endB):
    ''' 
        Minimal distance of two lines (here the vectors spanned by TMD ends)
        is calculated as the distance between parallel planes of the two lines
    '''
    ### Calculate determinant of vectors Aa-Ae=A and Ba-Be=B
    W1 = np.dot((begA - begB), (endA - begA))
    W2 = np.dot((begA - begB), (endB - begB))
    U11 = (np.linalg.norm(endA - begA)) ** 2
    U12 = np.dot((endA - begA), (endB - begB))
    U22 = (np.linalg.norm(endB - begB)) ** 2
    Det = np.dot(U11, U22) - np.dot(U12, U12) # The determinant is the area spanned by AxB
    
    if Det == 0:
        SA = 0.5
        SB = 0.5

    else:
        SA = (np.dot(W2, U12) - np.dot(W1, U22)) / Det
        SB = (np.dot(W2, U11) - np.dot(W1, U12)) / Det

        if SA > 1:
            SA = 1
        elif SA < 0:
            SA = 0
        else:
            SA = SA

        if SB > 1:
            SB = 1
        elif SB < 0:
            SB = 0
        else:
            SB = SB

    tA = begA + SA * (endA - begA)
    tB = begB + SB * (endB - begB)
    tA = np.array(tA, dtype=float)
    tB = np.array(tB, dtype=float)

    return tA, tB


