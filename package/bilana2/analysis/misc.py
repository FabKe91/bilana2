'''
    INCOMPLETE
'''

import logging
import MDAnalysis as mda

LOGGER = logging.getLogger("bilana2.analysis.misc")


def get_distance(systeminfo, reference_atoms: mda.AtomGroup, reference_identifier, selstr, selection_identifier, outputname=None, dim=2):
    ''' '''
    outformat_str = "{: <20}{: <15}{: <15}{: <20}{: <15}{: >15}{: >15}{: >15}"
    outformat = "{: <20}{: <15}{: <15}{: <20}{: <15}{: >15.5f}{: >15.5f}{: >15.5f}"
    u = systeminfo.universe
    len_traj = len(u.trajectory)
    sel     = u.select_atoms(selstr)
    if outputname is None:
        outputname = "distance_ref{}_sel{}_dim{}.dat".format(reference_identifier, selection_identifier, dim).replace(" ", "_")
    with open(outputname, "w") as outf:
        print(outformat_str.format("time", "resid", "resname", "atom", "distance", "x", "y", "z"), file=outf)
        for t in range(len_traj):
            time = u.trajectory[t].time
            if systeminfo.t_end < time or systeminfo.t_start > time:
                continue
            LOGGER.info("Time %s", time)

            pos_sel = sel.positions
            pos_ref = reference_atoms.center_of_mass()

            box     = u.dimensions
            if dim == 2:
                pos_ref[2] = 0
                pos_sel[:,2] = 0
                box[2] = 1
            elif dim == 3:
                pass
            else:
                raise ValueError("Dimension not possible. Use either 2 or 3")
            LOGGER.debug("Sel %s", sel,)
            LOGGER.debug("Resids %s", sel.resids)
            pos_array = mda.lib.distances.distance_array(pos_ref, pos_sel, box=box)[0]
            LOGGER.debug("Array: %s", pos_array)
            for ind, res in enumerate(sel.resids):
                resname = sel[ind].resname
                atom = sel[ind].name
                pos = sel[ind].position
                dist = pos_array[ind] / 10
                print(outformat.format(time, res, resname, atom, dist, *pos), file=outf)
