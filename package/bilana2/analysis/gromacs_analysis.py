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
    systeminfo.path.create_folders()
    os.makedirs(systeminfo.path.data + "/densities", exist_ok=True)
    TRJ = systeminfo.path.trj
    TPR = systeminfo.path.tpr
    NDX = systeminfo.path.tmp + "/temp.ndx"
    OUT = systeminfo.path.data + "/densities/" + outname
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
    commandstring = "-f {} -s {} -n {} -o {}  -d Z".format(TRJ, TPR, NDX, OUT)
    cmd = commandstring.split() + additional_input
    out, err = exec_gromacs(GMXNAME, "density", cmd)
    write_log("gmx_density", out, err, path=systeminfo.path.log)

    return OUT


def calc_rdf(systeminfo, ref, sel,
    seltype='atom', selrpos='atom', binsize=0.002, refprot=False, **kw_rdf):
    ''' Calculates 2D RDF of sel relative to ref only for one specific leaflet
        leaflet_assignment.dat (created in leaflets.py) is needed
        Using gmx rdf
    '''

    LOGGER.info("Calculating radial distribution function")
    LOGGER.info("Ref: %s\nSel: %s\n", ref, sel)
    systeminfo.path.create_folders()
    os.makedirs(systeminfo.path.data+'/rdf', exist_ok=True)

    ### Prepare addition keywords to fit exec_gromacs function ###
    additional_commands = []
    for key, val in kw_rdf.items():
        additional_commands += ["-"+key, str(val)]

    for leafndx, reslist in enumerate(systeminfo.convert.leaflet_resids):

        ### Create temporary files containing the selection string for -ref and -sel flag ###
        resid_list_str = 'resid ' + ' '.join([str(i) for i in reslist])
        selectdict = {}
        for selection in (('ref', ref), ('sel', sel)):

            if selection[0] == "ref" and refprot: # If protein is reference take all resids
                selectstring = '{}'.format(selection[1])
            else:
                selectstring = '{} {}'.format(selection[1], resid_list_str)

            select_fname = '{}/select{}_{}'.format(systeminfo.path.tmp,
                                                   selection[0], selection[1]).replace(" ", "_")
            selectdict.update({selection[1]:select_fname})

            with open(select_fname, "w") as selfile:
                print("Selection for {} is:\n{}\n".format(selection[0], selectstring))
                print(selectstring, file=selfile)

        ### Define file paths ###
        outputfile = '{}/rdf/rdf_leaf{}_{}-{}.xvg'.format(systeminfo.path.data, leafndx,
                                                   ref, sel).replace(" ", "")
        outputfile_cn = '{}/rdf/nr_leaf{}_{}-{}.xvg'.format(systeminfo.path.data, leafndx,
                                                     ref, sel).replace(" ", "")

        ### Run gromacs command ###
        cmd = [
            '-xy', '-xvg', 'none',
            '-f', systeminfo.path.trj, '-s', systeminfo.path.tpr,
            '-o', outputfile,  '-cn', outputfile_cn,
            '-ref',  '-sf', selectdict[ref],
            '-sel',  '-sf', selectdict[sel],
            '-selrpos', selrpos, '-seltype', seltype, '-bin', str(binsize),
            ]
        cmd += additional_commands
        out, err = exec_gromacs(GMXNAME, "rdf", cmd)
        write_log("gmx_rdf", out, err, path=systeminfo.path.log)
