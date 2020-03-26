import logging
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import distance_array

LOGGER = logging.getLogger("bilana2.analysis.interleaflet")

def calculate_overlapscore(sysinfo, refsel, chainsel,
    radius=10, delta_t=1000, nbins=25,
    outputfilename="overlap_score.csv",
    ):
    '''
        Calculate the overlap of chains within the bilayer.

        The scoring factor is based on the Bhattacharyya coeffectient
            ## BC(p,q) = norm * sum_i( sqrt(p_i*q_i) )
            ## norm = sqrt( 1/( sum(p) * sum(q) ) )

        refsel:     Selection of reference atoms to give a reference position
        chainsel:   Selection of atoms to calculate the overlap fore
        radius:     Gives the radius of the chunk size around each reference lipid
        delta_t:    Time frame to average over in ps
        nbins:      Number of bins to discretize z range of the bilayer

        Output layout of file is:
            <time> <resid> <thickness> <overlap> <chunksize>
    '''
    resids = []
    overlaps = []
    chunksizes = []
    chunkframes = []
    u = sysinfo.universe

    if len(u.residues.atoms.select_atoms(refsel)) == 0:
        raise ValueError('Refsel "{}" never matches any atoms'.format(refsel))
    if len(u.residues.atoms.select_atoms(chainsel)) == 0:
        raise ValueError('chainsel "{}" never matches any atoms'.format(refsel))

    for ts in u.trajectory:
        time = ts.time
        if not sysinfo.within_timerange(time):
            continue
        LOGGER.info("At %s ps", time)

        box      = u.dimensions
        residues = u.residues

        all_pos         = u.residues.atoms.select_atoms(refsel).positions
        all_pos_2d      = all_pos.copy()
        all_pos_2d[:,2] = 0

        for res_i in residues:
            if res_i.resname not in sysinfo.molecules or not (res_i.atoms.select_atoms(refsel)):
                continue

            ### Get positions for residue i ###
            ref_pos         = res_i.atoms.select_atoms(refsel).atoms.positions
            if ref_pos.shape[0] > 1:
                raise ValueError("only one reference atom per molecule allowed. Found more for"\
                    "res: {}".format(res_i))
            ref_pos = ref_pos[0]
            ref_pos_2d      = ref_pos.copy()
            ref_pos_2d[2]   = 0

            ### Define chunks per residue ###
            d_array          = distance_array(ref_pos_2d, all_pos_2d, box=box)[0]
            within_chunk     = np.where(d_array <= radius)
            res_within_chunk = residues[within_chunk]
            chunksize        = len(res_within_chunk)

            chunk_pos = res_within_chunk.atoms.select_atoms(refsel).positions
            chunk_pos_z_mean = chunk_pos[:,2].mean()

            res_lower = res_within_chunk.atoms.select_atoms(
                    refsel+" and prop z <= {}".format(chunk_pos_z_mean) ).residues
            res_upper = res_within_chunk.atoms.select_atoms(
                    refsel+" and prop z >  {}".format(chunk_pos_z_mean) ).residues

            ### Calculate Score value ####
            tail_pos       = res_within_chunk.atoms.select_atoms(chainsel).positions
            min_z = min(tail_pos[:,2])
            max_z = max(tail_pos[:,2])
            tail_pos_lower = res_lower.atoms.select_atoms(chainsel).positions
            tail_pos_upper = res_upper.atoms.select_atoms(chainsel).positions

            zdistr_upper = np.histogram( tail_pos_upper[:,2],
                        range=(min_z, max_z), bins=nbins )[0]
            zdistr_lower = np.histogram( tail_pos_lower[:,2],
                        range=(min_z, max_z), bins=nbins )[0]

            sum_distr = np.sqrt( np.sum(zdistr_upper) * np.sum(zdistr_lower) )

            score = sum_distr**-1 * np.sum( np.sqrt( zdistr_upper * zdistr_lower ) )

            ### Add data to container ###
            resids.append(res_i.resid)
            overlaps.append(score)
            chunksizes.append(chunksize)

        chunkframe = pd.DataFrame({"resid":resids, "overlap":overlaps, "chunksize":chunksizes})
        chunkframe["time"] = time
        chunkframes.append( chunkframe )

        LOGGER.debug("from residue choice:\n %s \n %s \n", res_upper, res_lower)
        LOGGER.debug("Last entries of chunk:")
        LOGGER.debug("res %s", resids[-1])
        LOGGER.debug("overlap %s, with distr %s", overlaps[-1], zdistr_lower)
        LOGGER.debug("chunksize %s", chunksizes[-1])

    final = pd.concat(chunkframes)

    LOGGER.debug("Frame before averaging:\n%s", final)

    ### Calculating number of frames to be averaged over ###
    dt         = u.trajectory.dt
    frame_step = int(delta_t // dt) if int(delta_t // dt) else 1

    if frame_step != 1: ### No need of averaing if delta_t already time step of simulation
        start_time = frame_step * dt
        final_time = len(u.trajectory) - frame_step * dt
        binrange = (start_time, final_time, delta_t)
        bins =  []
        for i in np.arange(*binrange):
            bins.append( (i,i+binrange[2]) )
        bins = pd.IntervalIndex.from_tuples(bins)
        final["bins"] = pd.cut(final.time, bins=bins)
        final = final.groupby(["bins", "resid"]).mean().reset_index().drop(columns=["bins"])

    final.to_csv(outputfilename, index=False)
    return outputfilename

def calculate_thickness(sysinfo, refsel,
    radius=10, delta_t=1000,
    outputfilename="thickness_chunks.csv",
    ):
    '''
        Calculate the thickness of chunks in a bilayer.
        A chunk is comprised of all lipids within a cylinder around each residue with radius r
        (using only xy distance)

        refsel:     Selection of reference atoms to give a reference position
        radius:     Gives the radius of the chunk size around each reference lipid
        delta_t:    Time frame to average over in ps

        Output layout of file is:
            <time> <resid> <thickness> <chunksize>
    '''
    resids = []
    chunksizes = []
    thicknesses = []
    chunkframes = []
    u = sysinfo.universe

    if not len(u.residues.atoms.select_atoms(refsel)) == 0:
        raise ValueError('refsel "{}" never matches any atoms'.format(refsel))

    for ts in u.trajectory:
        time = ts.time
        LOGGER.info("At %s ps", time)
        if not sysinfo.within_timerange(time):
            continue

        ### Doing the calculation for each chunk ###
        box      = u.dimensions
        residues = u.residues

        all_pos         = u.residues.atoms.select_atoms(refsel).positions
        all_pos_2d      = all_pos.copy()
        all_pos_2d[:,2] = 0

        for res_i in residues:
            if res_i.resname not in sysinfo.molecules or not (res_i.atoms.select_atoms(refsel)):
                continue

            ### Get positions for residue i ###
            ref_pos         = res_i.atoms.select_atoms(refsel).atoms.positions
            if ref_pos.shape[0] > 1:
                raise ValueError("only one reference atom per molecule allowed. Found more for"\
                    "res: {}".format(res_i))
            ref_pos = ref_pos[0]
            ref_pos_2d      = ref_pos.copy()
            ref_pos_2d[2]   = 0

            ### Define chunks per residue ###
            d_array          = distance_array(ref_pos_2d, all_pos_2d, box=box)[0]
            within_chunk     = np.where(d_array <= radius)
            res_within_chunk = residues[within_chunk]
            chunksize        = len(res_within_chunk)

            chunk_pos = res_within_chunk.atoms.select_atoms(refsel).positions

            chunk_pos_z_mean = chunk_pos[:,2].mean()

            chunk_pos_lower = chunk_pos[chunk_pos[:,2] <= chunk_pos_z_mean]
            chunk_pos_upper = chunk_pos[chunk_pos[:,2]  > chunk_pos_z_mean]

            ### Now calculating thickness ###
            chunk_pos_mean_upper = chunk_pos_upper[:,2].mean()
            chunk_pos_mean_lower = chunk_pos_lower[:,2].mean()

            thickness = np.abs( chunk_pos_mean_upper - chunk_pos_mean_lower )

            ### Add data to container ###
            resids.append(res_i.resid)
            thicknesses.append(thickness)
            chunksizes.append(chunksize)

        chunkframe = pd.DataFrame({"resid":resids, "thickness":thicknesses, "chunksize":chunksizes})
        chunkframe["time"] = time
        chunkframes.append( chunkframe )

        LOGGER.debug("from residue choice:\n %s \n", res_within_chunk)
        LOGGER.debug("Last entries of chunk:")
        LOGGER.debug("res %s", resids[-1])
        LOGGER.debug("thickness %s, from pos: %s - %s",
            thicknesses[-1], chunk_pos_mean_upper, chunk_pos_mean_lower )
        LOGGER.debug("chunksize %s", chunksizes[-1])


    final = pd.concat(chunkframes)

    LOGGER.debug("frame before averaging:\n%s", final)

    ### Calculating number of frames to be averaged over ###
    dt         = u.trajectory.dt
    frame_step = int(delta_t // dt) if int(delta_t // dt) else 1

    if frame_step != 1:
        start_time = frame_step * dt
        final_time = len(u.trajectory) - frame_step * dt
        binrange = (start_time, final_time, delta_t)
        bins =  []
        for i in np.arange(*binrange):
            bins.append( (i,i+binrange[2]) )
        bins = pd.IntervalIndex.from_tuples(bins)
        final["bins"] = pd.cut(final.time, bins=bins)
        final = final.groupby(["bins", "resid"]).mean().reset_index().drop(columns=["bins"])

    final.to_csv(outputfilename, index=False)
    return outputfilename