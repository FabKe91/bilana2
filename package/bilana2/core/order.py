'''
    ================================================================
    Core objects --- :mod: 'bilana.core.order :
        create_cc_orderfiles
        calc_tilt
    ================================================================


'''
import logging
import numpy as np
import pandas as pd
from ..lib.common import angle_to_axis

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

LOGGER = logging.getLogger("bilana2.core.order")

def create_cc_orderfiles(sysinfo,
    outputfile_scd='scd_distribution.dat',
    outputfile_s_profile='s_profile.dat',
    with_tilt_correction='tilt.csv'):
    '''

        Two files will be written:
            - outputfile_scd:       Contains the average order parameter value per lipid
                                    at each time
            - outputfile_s_profile: Stores the order parameter profiles per lipid molecule
                                    at each time

        The order parameter S = 0.5 (3cos^2(X)-1) is calculated in the function
        get_cc_order()

        X is the angle between the bilayer normal and (set with mode):
        the vector Cn-Cn+2 carbon atoms of chain. Averaged over all the angles

        Tilt correction will read the average tilt vector per leaflet at time and use it as
        reference axis

    '''
    # If tilt correction is activated read file with tilt information or create it
    if with_tilt_correction:
        try:
            dat = pd.read_csv(with_tilt_correction)
        except FileNotFoundError:
            LOGGER.info("Could not find tilt file. Creating new one.")
            calc_tilt(sysinfo)
            dat = pd.read_csv(with_tilt_correction)

        # Define how the order parameter has to be calculated

    ## Gather all input data for _calc_scd_output function
    len_traj = len(sysinfo.universe.trajectory)

    with open(outputfile_scd, "w") as scdfile, open(outputfile_s_profile, "w") as sprof_file:

        #### Print header files ####
        print("{: <12}{: <10}{: <10}{: <7}{: <15}"\
                .format("time", "resid", "leaflet", "resname", "Scd"),
            file=scdfile)
        print("{: <12}{: <10}{: <10}{: <7}{: <15}{: <15}{: <10}{: <10}"\
            .format("time", "resid", "leaflet", "resname", "avgS", "S", "carbon", "chain" ),
            file=sprof_file
            )

        for t in range(len_traj):

            time = sysinfo.universe.trajectory[t].time
            if not sysinfo.within_timerange(time):
                continue
            LOGGER.info("At time %s", time)

            ### Get tilt of leaflet at time ###
            if with_tilt_correction:
                new_axis = np.asarray(dat.loc[(dat.time == time)][["x", "y", "z"]]).copy()
                LOGGER.debug("Corrected angle %s", new_axis)
            else:
                new_axis = None

            for res in sysinfo.lipid_resids:
                LOGGER.debug("At time %s and residue %s", time, res)

                leaflet = sysinfo.convert.resid_to_leaflet[res]
                resname = sysinfo.convert.resid_to_resname[res]

                if new_axis is not None:
                    new_axis_at_t = new_axis[leaflet]

                if resname not in sysinfo.molecules:
                    continue

                LOGGER.debug("getting to positions ...")
                ### Get positions ###
                tailatms = sysinfo.ff.scd_tail_atoms_of(resname)
                positions = []
                for sn in tailatms:
                    pos = sysinfo.universe.atoms.select_atoms( "resid {} and name {}"\
                            .format(res, ' '.join(sn) ) ).positions
                    positions.append(pos)

                order_val, s_prof = get_cc_order(positions, ref_axis=new_axis_at_t)

                LOGGER.debug("printing to files ...")
                ### Print everything to files ###
                line_scd = "{: <12.2f}{: <10}{: <10}{: <7}{: <15.8}".format(
                        time, res, leaflet, resname, order_val)
                print(line_scd, file=scdfile)
                for chain_ndx, slist in enumerate(s_prof):
                    for carb_i, order_carb in enumerate(slist):
                        line_p = "{: <12.2f}{: <10}{: <10}{: <7}{: <15.8}{: <15.8}{: <10}{: <10}"\
                            .format(time, res, leaflet, resname,
                            order_val, order_carb, carb_i, chain_ndx)
                    print(line_p, file=sprof_file)

def get_cc_order(positions: [np.array,], ref_axis=(0,0,1)) -> float:
    ''' Calculate the cc order parameter

        cc order parameter is defined as 0.5 * ( ( 3 * (cos_angle**2)) - 1 )
        for each consecutive positions in the positions array, whole the cos_angle is
        calculated with respect to ref_axis

        Input positions must be a list of arrays of positions:
        The list holds the arrays containing atom positions
        while the first array holds positions of chain 1, the second of chain 2 ...
        e.g.: np.array( np.array(pos_sn1), np.array(pos_sn2) )


    '''
    assert isinstance(positions, list)

    s_vals = [[] for _ in positions]

    for sn_x, positions_sn_x in enumerate(positions):

        #scds_of_atoms = []
        for i in range( len(positions_sn_x) - 1 ):# Explicitly using range(len()) to save if clause

            pos1, pos2 = positions_sn_x[i], positions_sn_x[i+1]

            diffvector = pos2 - pos1
            diffvector /= np.linalg.norm(diffvector)

            cos_angle = np.dot(diffvector, ref_axis)
            #scds_of_atoms.append( 0.5 * ( ( 3 * (cos_angle**2)) - 1 )  )
            s_vals[sn_x].append( 0.5 * ( ( 3 * (cos_angle**2)) - 1 )  )

            LOGGER.debug("Diffvector %s", diffvector)
            LOGGER.debug("Resulting cos %s", cos_angle)
        s_vals[sn_x] = np.array( s_vals[sn_x] )

    return np.array(s_vals).mean(), s_vals


def calc_tilt(sysinfo, filename="tilt.csv"):
    '''
        Calculate the tilt angle per _PL_ molecule in each leaflet with regard to the z axis
    '''
    u = sysinfo.universe
    len_traj = len(u.trajectory)

    angles              = []
    vectors             = []
    times               = []
    leaflet_assignments = []

    for t in range(len_traj):

        time = u.trajectory[t].time
        if not sysinfo.within_timerange(time):
            continue

        LOGGER.info("At time %s", time)

        leaflist = [[], []]

        for res in sysinfo.lipid_resids:
            LOGGER.debug("at res %s", res)

            resname = sysinfo.convert.resid_to_resname[res]
            leaflet = sysinfo.convert.resid_to_leaflet[res]

            if resname not in sysinfo.pl_resnames:
                continue

            ### first carbon of tails ###
            masks1 = ( u.atoms.names == sysinfo.ff.tailcarbons_of(resname)[0][0],
                    u.atoms.resids == res )
            masks2 = ( u.atoms.names == sysinfo.ff.tailcarbons_of(resname)[1][0],
                    u.atoms.resids == res )

            t1_xyz = u.atoms.positions[ masks1[0] & masks1[1] ][0]
            t2_xyz = u.atoms.positions[ masks2[0] & masks2[1] ][0]

            ### get vector c_last - c_first ###
            masks1 = ( u.atoms.names == sysinfo.ff.tailcarbons_of(resname)[0][-1],
                    u.atoms.resids == res )
            masks2 = ( u.atoms.names == sysinfo.ff.tailcarbons_of(resname)[1][-1],
                    u.atoms.resids == res )

            tail1_xyz = u.atoms.positions[ masks1[0] & masks1[1] ][0] - t1_xyz
            tail2_xyz = u.atoms.positions[ masks2[0] & masks2[1] ][0] - t2_xyz
            ### normalize vectors ###
            tail1_xyz /= np.linalg.norm(tail1_xyz)
            tail2_xyz /= np.linalg.norm(tail2_xyz)

            LOGGER.debug("tail1_xyz %s", tail1_xyz)

            ### leaf is either 0 or 1: append avg vector to correct leaflet list ###
            leaflist[leaflet] += [tail1_xyz, tail2_xyz]
        LOGGER.debug("leaflist: %s", leaflist)

        ### get average tilt vector per leaflet ###
        avg_vec_leaf1 = np.array(leaflist[0]).mean(axis=0)
        avg_vec_leaf2 = np.array(leaflist[1]).mean(axis=0)
        avg_vec_leaf1 = (avg_vec_leaf1/np.linalg.norm(avg_vec_leaf1))
        avg_vec_leaf2 = (avg_vec_leaf2/np.linalg.norm(avg_vec_leaf2))
        LOGGER.debug("avg_vec_leaf1 %s", avg_vec_leaf1)

        ### calculate the angle of average vector to z axis ###
        ang_l1 = angle_to_axis(avg_vec_leaf1)
        LOGGER.debug("ang1 %s", ang_l1)
        ang_l2 = angle_to_axis(avg_vec_leaf2)
        if ang_l1 > 90:
            ang_l1 = np.abs(ang_l1 - 180)
        if ang_l2 > 90:
            ang_l2 = np.abs(ang_l2 - 180)

        ### fill lists this will become lines in the dataframe ###
        angles              += [ang_l1, ang_l2]
        vectors             += [avg_vec_leaf1, avg_vec_leaf2]
        times               += [time, time]
        leaflet_assignments += [0, 1]

    vectors = np.array(vectors)

    dat = pd.DataFrame({"time":times, "leaflet":leaflet_assignments, "angle":angles,
            "x":vectors[:,0], "y":vectors[:,1], "z":vectors[:,2]})
    dat = dat.sort_values(by=["time"])

    dat.to_csv(filename, index=False)