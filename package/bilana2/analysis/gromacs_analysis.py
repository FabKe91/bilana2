import os
import logging

from ..lib.gromacswrapper import exec_gromacs, GMXNAME, write_log

LOGGER = logging.getLogger("bilana2.analysis.gromacs_analysis.py")

def calc_density(systeminfo, selstr, outname="density.xvg", overwrite=False, **kw_den):
    '''
        Uses density calculation of gromacs
        1. Get index file using gmx select with
            gmx select -f ... -select <selstr>
        2. Run
            gmx density -f ... -d Z
            NOTE: Additional flags can be set adding with kw_den like b=3 converted to -b 3
    '''
    ### Setting all paths needed ####
    os.makedirs(systeminfo.path.data + "/densities", exist_ok=True)
    TRJ = systeminfo.path.trj
    TPR = systeminfo.path.tpr
    NDX = systeminfo.path.tmp + "/temp.ndx"
    OUT = systeminfo.path.datapath + "/densities/" + outname
    if not overwrite and os.path.exists(OUT):
        LOGGER.warning("Density file already exists")
        return OUT

    ### Preparing additional input to conform exec_gromacs function ###
    additional_input = []
    if kw_den:
        keys = kw_den.keys()
        keys = ["-"+i for i in keys]
        vals = kw_den.values()
        for z in zip(keys, vals):
            additional_input += list(z)

    ####
    LOGGER.info("Creating index file...")
    ####
    commandstring = '-f {} -s {} -on {} -select'.format(TRJ, TPR, NDX) ## dont forget the selstr
    cmd = commandstring.split() + [selstr]
    out, err = exec_gromacs(GMXNAME, "select", cmd)
    write_log("gmx_select", out, err, path=systeminfo.path.log)

    ####
    LOGGER.info("Run gmx density...")
    ####
    commandstring = "gmx density -f {} -s {} -n {} -o {}  -d Z".format(TRJ, TPR, NDX, OUT)
    cmd = commandstring.split() + additional_input
    out, err = exec_gromacs(GMXNAME, "density", cmd)
    write_log("gmx_density", out, err, path=systeminfo.path.log)

    return OUT