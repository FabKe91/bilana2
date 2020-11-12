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
    commandstring = "-f {} -s {} -n {} -o {}  -d Z -b {} -e {}".format(TRJ, TPR, NDX, OUT, systeminfo.t_start, systeminfo.t_end)
    cmd = commandstring.split() + additional_input
    out, err = exec_gromacs(GMXNAME, "density", cmd)
    write_log("gmx_density", out, err, path=systeminfo.path.log)

    return OUT


def calc_rdf(systeminfo, ref, sel,
    seltype='atom', selrpos='atom', 
    binsize=0.002, protein=None, name_suffix="",
    **kw_rdf):
    ''' Calculates 2D RDF of sel relative to ref only for one specific leaflet
        leaflet_assignment.dat (created in leaflets.py) is needed
        Using gmx rdf
        
        protein must be Protein() class instance
        if specified: will be used as reference
            only those AA in same leaflet will be used 

    '''
    if name_suffix and name_suffix[0] != "_":
        name_suffix = "_"+name_suffix
    if protein is not None:
        if isinstance(protein, list):
            protidname = "_prot{}".format('-'.join([str(prot.id) for prot in protein]))
        else:
            protidname = "_prot{}".format(protein.id)
    else: 
        protidname = ""

    LOGGER.info("Calculating radial distribution function")
    LOGGER.info("Ref: %s\nSel: %s\n", ref, sel)
    systeminfo.path.create_folders()
    os.makedirs(systeminfo.path.data+'/rdf', exist_ok=True)

    ### Prepare addition keywords to fit exec_gromacs function ###
    additional_commands = []
    for key, val in kw_rdf.items():
        additional_commands += ["-"+key, str(val)]

    for leafndx, reslist in enumerate(systeminfo.convert.leaflet_to_resids):

        ### Create temporary fies containing the selection string for -ref and -sel flag ###
        resid_list_str = 'resid ' + ' '.join([str(i) for i in reslist])
        
        pl_in_leaf = True

        selectdict = {}
        for selection in (('ref', ref), ('sel', sel)):

            # if protein is used as ref: use only residues in similar leaflet as leafndx
            if protein is not None: 

                #if both ref and sel are protein
                if isinstance(protein, list):

                    if selection[0] == "ref":
                        protndx = 0
                    elif selection[0] == "sel":
                        protndx = 1

                    prot_residlist_str = 'resid ' + ' '.join(
                            [ str(protres) for protres, protleafndx in protein[protndx].resid_to_leaflet.items()\
                                if protleafndx == leafndx ] )

                    LOGGER.debug("Protein residlist %s", prot_residlist_str)
                    selectstring = '{} and {}'.format(selection[1], prot_residlist_str)

                ### if only ref is protein
                else:
                    if selection[0] == "ref":

                        prot_residlist_str = 'resid ' + ' '.join(
                                [ str(protres) for protres, protleafndx in protein.resid_to_leaflet.items()\
                                    if protleafndx == leafndx ] )

                        LOGGER.debug("Protein residlist %s", prot_residlist_str)
                        selectstring = '{} and {}'.format(selection[1], prot_residlist_str)
                    elif selection[0] == "sel":
                        selectstring = '{} and {}'.format(selection[1], resid_list_str)
                        print([lip for lip in systeminfo.universe.atoms.select_atoms(selection[1]) if systeminfo.convert.resid_to_leaflet[lip.residue.resid] == leafndx ])

                        if not [lip for lip in systeminfo.universe.atoms.select_atoms(selection[1]) if systeminfo.convert.resid_to_leaflet[lip.residue.resid] == leafndx ]:
                            pl_in_leaf = False

            ### if both ref and sel are lipids
            else:
                selectstring = '{} and {}'.format(selection[1], resid_list_str)
                print([lip for lip in systeminfo.universe.atoms.select_atoms(selection[1]) if systeminfo.convert.resid_to_leaflet[lip.residue.resid] == leafndx ])

                if not [lip for lip in systeminfo.universe.atoms.select_atoms(selection[1]) if systeminfo.convert.resid_to_leaflet[lip.residue.resid] == leafndx ]:
                    pl_in_leaf = False
                
            select_fname = '{}/select{}_{}'.format(systeminfo.path.tmp,
                                                   selection[0], selection[1]).replace(" ", "_")
            selectdict.update({selection[1]:select_fname})

            with open(select_fname, "w") as selfile:
                print("Selection for {} is:\n{}\n".format(selection[0], selectstring))
                print(selectstring, file=selfile)

        if not pl_in_leaf:
            print("skipping", selection)
            continue

        ### Define file paths ###
        outputfile = '{}/rdf/rdf_leaf{}{}_{}-{}{}.xvg'.format(systeminfo.path.data,
                                                   leafndx, protidname,
                                                   ref, sel, name_suffix)\
                                                   .replace(" ", "").replace("\"", "")
        outputfile_cn = '{}/rdf/nr_leaf{}{}_{}-{}{}.xvg'.format(systeminfo.path.data, 
                                                   leafndx, protidname,
                                                   ref, sel, name_suffix)\
                                                   .replace(" ", "").replace("\"", "")

        ### Run gromacs command ###
        cmd = [
            '-xy', '-xvg', 'none',
            '-f', systeminfo.path.trj, '-s', systeminfo.path.tpr,
            '-o', outputfile,  '-cn', outputfile_cn,
            '-ref',  '-sf', selectdict[ref],
            '-sel',  '-sf', selectdict[sel],
            '-selrpos', selrpos, '-seltype', seltype, '-bin', str(binsize),
            '-b', str(systeminfo.t_start), '-e', str(systeminfo.t_end),
            ]
        cmd += additional_commands
        out, err = exec_gromacs(GMXNAME, "rdf", cmd)
        write_log("gmx_rdf", out, err, path=systeminfo.path.log)


def calc_diffusion(sysinfo, selection, mol=False):
    ''' Use gromacs to calculate diffusion of lipid type '''
    outpath = sysinfo.path.data + "/diffusion/"
    os.makedirs(outpath, exist_ok=True)

    # get index entry for atom choice
    ndxname = sysinfo.path.tmp + "/tmp{}.ndx".format(selection.replace(" ", ""))
    selstr = 'sel={};sel'.format(selection)
    cmd = [ 
        "-s", sysinfo.path.gro, "-on", ndxname,
        "-select", selstr,
    ]
    out, err = exec_gromacs(GMXNAME, "select", cmd)
    write_log("gmx_select", out, err, path=sysinfo.path.log)

    # Calculate MSD and diffusion
    output_msd = outpath + "/msd_{}.xvg".format(selection.replace(" ", "").replace("\"", ""))
    output_diff = outpath + "/diff_{}.xvg".format(selection.replace(" ", "").replace("\"", ""))
    inpstr = "sel" + "\n"
    cmd = [
        "-f", sysinfo.path.trj, "-s", sysinfo.path.tpr, "-n", ndxname,
        "-o", output_msd,
        "-lateral", "z", "-xvg", "none"
    ]
    if mol:
        cmd += ["-mol", output_diff]
    out, err = exec_gromacs(GMXNAME, "msd", cmd, inpstr)
    write_log("gmx_msd", out, err, path=sysinfo.path.log)
    

