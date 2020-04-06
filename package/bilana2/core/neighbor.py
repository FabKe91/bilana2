'''
    ==========================================================================================
    Core objects: --- :mod: 'bilana.core.neighbor
        - write_neighbor_info   Creates a file that searches all neighbors of each
                                molecule per frame
        - get_neighbor_dict     Reads neighbor file created in write_neighbor_info and returns
                                a dictionary for easy handling
    ==========================================================================================
'''

import logging
import numpy as np
import pandas as pd
import MDAnalysis as mda
#from ..lib.common import molecule_leaflet_orientation

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

LOGGER = logging.getLogger("bilana2.core.neighbor")

def write_neighbor_info(sysinfo, outputfilename="neighbor_info", mode="atom"):
    ''' Creates neighbor_info file based on 2d distances (xy plane) '''

    refatoms         = sysinfo.reference_atomselection

    with open(outputfilename, "w") as outf:
        print("{: <20}{: <20}{: <20}{: <20}"\
            .format("Resid", "Time", "Number_of_neighbors", "List_of_Neighbors"), file=outf)

        traj_len = len(sysinfo.universe.trajectory)
        for t in range(traj_len):

            time = sysinfo.universe.trajectory[t].time
            if not sysinfo.within_timerange(time):
                continue
            LOGGER.info("At time %s", time)

            refatomgrp = sysinfo.universe.select_atoms(refatoms)
            refpositions = get_ref_positions(sysinfo, mode, refatomgrp) # leaflets=[(resid1, pos1), ...]

            LOGGER.debug("Found %s atoms: %s ", len(refatomgrp.atoms), refatomgrp.atoms)
            LOGGER.debug("Dimension of leaflets %s", np.array(refpositions).shape)

            ### Leaflet thickness is the standarddeviation of z positions ###
            ### of heads around bilayer center ###
            leaflet_thickness = np.sqrt(refatomgrp.atoms.positions[:,2].var())

            for hostid, host_pos in refpositions:

                host_pos2d = host_pos.copy()
                host_pos2d[2] = 0

                position_array   =  np.array([pos for resid, pos in refpositions])
                position_array2d =  position_array.copy()
                position_array2d[:,2] = 0

                dist_array = mda.lib.distances.distance_array(
                    host_pos,
                    position_array,
                    box=sysinfo.universe.dimensions)[0]
                dist_array2d = mda.lib.distances.distance_array(
                    host_pos2d,
                    position_array2d,
                    box=sysinfo.universe.dimensions)[0]

                neiblist = []
                for resndx, distance in enumerate(dist_array):

                    if distance <= leaflet_thickness:
                        distance2d = dist_array2d[resndx]
                        if distance2d <= sysinfo.cutoff*10.0:
                            neiblist.append( refpositions[resndx][0] )

                neiblist = list(set(neiblist)) # delete duplicates
                neiblist.sort()
                neiblist = [ str(n) for n in neiblist if n != hostid ] # delete host entry
                n_neibs = len(neiblist)

                line = "{: <20}{: <20}{: <20}{: <20}"\
                        .format( hostid, time, n_neibs, ','.join(str(i) for i in neiblist) )
                print(line, file=outf)


def get_ref_positions(sysinfo, mode, refatomgrp):
    '''
        Return the reference positions per leaflet and molecule depending on mode (see below)
        Possible modes
            - atom:     Use one atom per molecule as reference position and return a list of
                        two lists (per leaflet) with positions of all lipid molecules

            - center    !!Not yet implemented!! Use the center of geometry of an atom group
                        specified

            - tails     Use one atom per leaflet tail as reference position

        One critical aspect is the choice of a correct leaflet assignment - Right now the assignment
        is made from the first frame, while uncritical for phospholipids, this may be wrong due to
        sterol flip flop motion.
        In future sterol leaflet assignment should be dynamically
    '''

    modes = ["atom",]

    if mode == "atom":
        if len(refatomgrp.resids) != len(set(refatomgrp.resids)):
            raise ValueError("Refatoms string leads to more than one entry per molecule: {}-{}"\
                .format(refatomgrp, refatomgrp.resids))

        ### Get leaflet assignment dict for PLs and sterols ###
        #orientations = {}
        #for residue in refatomgrp.residues:
        #    #headsel = 'resname {} and name {}'\
        #       .format(residue.resname, ' '.join( .head_atoms_of(residue.resname) ) )
        #    #tailsel = 'resname {} and name {}'\
        #       .format(residue.resname, ' '.join( [taillist[-1] for taillist in\
        #           sysinfo.ff.tailcarbons_of(residue.resname) ] ) )
        #    #head_pos = residue.atoms.select_atoms( headsel ).center_of_mass()
        #    #tail_pos = residue.atoms.select_atoms( tailsel ).center_of_mass()
        #    if residue.resname in sysinfo.pl_resnames:
        #        orientations[residue.resid] = sysinfo.convert.resid_to_leaflet(residue.resid)
        #    elif:
        #        pass

        ### Get the actual reference positions per leaflet ###
        #leaf1 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms\
        #        if sysinfo.convert.resid_to_leaflet[ atm.resid ] ] )
        #leaf2 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms\
        #        if not sysinfo.convert.resid_to_leaflet[ atm.resid ] ] )
        refpositions = np.array( [ (atm.resid, atm.position,) for atm  in refatomgrp.atoms ] )

    #elif mode == "tails":
    #    cnt = 0
    #    orientations = {}
    #    for residue in refatomgrp.residues:
    #        headsel = 'resname {} and name {}'\
    #            .format(residue.resname, ' '.join( sysinfo.ff.head_atoms_of(residue.resname) ) )
    #        tailsel = 'resname {} and name {}'\
    #            .format(residue.resname,
    #                    ' '.join( np.array( sysinfo.ff.tailcarbons_of(residue.resname) )[:,-1] ) )
    #        LOGGER.debug("Head selection: %s", headsel)
    #        LOGGER.debug("Tail selection: %s", tailsel)
    #        head_pos = residue.atoms.select_atoms( headsel ).center_of_mass()
    #        tail_pos = residue.atoms.select_atoms( tailsel ).center_of_mass()
    #        orientations[residue.resid] =  molecule_leaflet_orientation( head_pos, tail_pos )
    #    leaf1 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms\
    #            if orientations[ atm.resid ] ] )
    #    leaf2 = np.array( [(atm.resid, atm.position) for atm  in refatomgrp.atoms\
    #            if not orientations[ atm.resid ] ] )
    #    LOGGER.debug("Leaf1 is\n%s", leaf1)
    #    for resid  in refatomgrp.resids:
    #        resid += cnt
    #        mask1 = np.where( leaf1[:,0] == resid )[0] # Index of duplicate per resid
    #        mask2 = np.where( leaf2[:,0] == resid )[0] # are stored here
    #        if mask1.shape[0] == 2:
    #            LOGGER.debug("Index of second tail: %s", np.where( leaf1[:,0] == resid ))
    #            cnt += 1 # count number of additional resids
    #            leaf1[ leaf1[:,0] > resid ] += 1 # increase resid number above processed resid
    #            leaf2[ leaf2[:,0] > resid ] += 1 # for each duplicate
    #            LOGGER.debug("leaf1 now %s", leaf1)
    #            leaf1[ mask1[1] ][0] += 1  # Increase resid of second item in resid duplicate
    #            LOGGER.debug("leaf1 and then %s", leaf1)
    #        elif mask2.shape[0] == 2: # Same here as for mask1
    #            LOGGER.debug("Index of second tail: %s", np.where( leaf2[:,0] == resid ))
    #            cnt += 1
    #            leaf1[ leaf1[:,0] > resid ] += 1
    #            leaf2[ leaf2[:,0] > resid ] += 1
    #            LOGGER.debug("leaf2 now %s", leaf2)
    #            leaf2[ mask2[1] ][0] += 1
    #            LOGGER.debug("leaf2 and then %s", leaf2)
    #    leaf1, leaf2 = list(leaf1), list(leaf2)

    else:
        raise ValueError("Invalid mode, choose one of {}".format(modes))
    return refpositions

def get_neighbor_dict(neighborfilename='neighbor_info'):
    ''' Returns a list of all neighbors being in the
        cutoff distance at least once in the trajectory.
        Neighborfile is required and is output of determine_neighbors function

        Dict layout is:
        neibdict[time][resid] -> [neibs]

    '''

    neibdict = {}

    data = pd.read_table(neighborfilename, delim_whitespace=True)
    data["nlist"] = data.fillna('').List_of_Neighbors\
        .apply(lambda x: [int(i) for i in x.split(',') if i ])
    data = data.drop( columns=["Number_of_neighbors", "List_of_Neighbors"] )

    for t, fr in data.groupby("Time"):
        neibdict[t] = fr.set_index("Resid").to_dict()["nlist"]

    return neibdict
