import logging
import numpy as np
import pandas as pd

LOGGER = logging.getLogger("bilana2.lib.files")

def create_eofs(outputfile="EofScd{}.csv",
    efile="all_energies.dat",
    sfile="scd_distribution.dat",
    neighborinfofile="neighbor_info",
    neighbortypefile="neighborcount.dat"):
    ''' Merges all_energies and scd_distribution file with output columns:
        var outputfile must contain a format control symbol ({}), so the interaction between
        specific lipid types is separated.

    <time> <host> <host_Scd> <neib> <neib_Scd> <molparts> <DeltaScd> <AvgScd> <leaflet> \
         <host_type> <neib_type> <resname1> <resname2> ... <resnameX>

    '''
    def compare_cols(expectedcol, col):
        if set(col) != set(expectedcol):
            raise ValueError("Columns do not fit!\nexpected: {}\nfound:{}".format(expectedcol, col))

    LOGGER.info("loading energy file (%s)...", efile)
    #HEAD: time   host   neighbor  molparts     VdW       Coul       Etot
    energy = pd.read_table(efile, sep="\s+")
    colnames = ["time", "host", "neighbor", "molparts", "VdW", "Coul", "Etot"]
    compare_cols(colnames, energy.columns)

    LOGGER.info("loading order file (%s)...", sfile)
    #HEAD: time   resid   leaflet  resname   Scd
    order  = pd.read_table(sfile, sep="\s+")
    colnames = ["time", "resid", "leaflet", "resname", "Scd"]
    compare_cols(colnames, order.columns)

    LOGGER.info("loading neighbortype file (%s)...", neighbortypefile)
    #HEAD: time   resid resname   Neibs_Type1 Neibs_Type2 ...
    neib = pd.read_table(neighbortypefile, sep="\s+")
    if not ("time" in neib.columns and "resid" in neib.columns and "resname" in neib.columns):
        compare_cols(["time", "resid", "resname", "..."], neib.columns)

    LOGGER.info("loading neighbor info file (%s)...", neighborinfofile)
    #HEAD: resid  time  Number_of_neighbors List_of_Neighbors
    neibmap = pd.read_table(neighborinfofile, sep="\s+")
    colnames = ["time", "resid", "Number_of_neighbors", "List_of_Neighbors"]
    compare_cols(colnames, neibmap.columns)


    LOGGER.info("editing names...")
    ### Edit names | From here on work with new column list ! ###
    neibmap["nlist"] = neibmap.fillna('').List_of_Neighbors\
            .apply(lambda x: [int(i) for i in x.split(',') if i ])
    neibmap = neibmap.drop( columns=["Number_of_neighbors", "List_of_Neighbors"] )
    neibmap = neibmap.rename(columns={"resid":"host"})
    neib   = neib.rename(columns={"resid":"host", "resname":"host_type"})
    order  = order.rename(columns={"resid":"host", "resname":"host_type", "Scd":"host_Scd"})
    energy = energy.rename(columns={"neighbor":"neib"})
    ### New column names are:
    ### energy:         time   host             neib      molparts     VdW       Coul       Etot
    ### order:          time   host   resname   leaflet  Scd
    ### neighbortype:   time   host   resname   Neibs_Type1 Neibs_Type2 ...
    ### neighbor info:  time   host             nlist

    ### !!! CHANGED THAT HERE !!! ###
    ### !!! NOW DUPLICATE ENTRIES ARE USED SO E_12 != E_21 !!! ###
    ####remove duplicate entries (E_12 == E_21) ###
    #LOGGER.info("remove duplicate entries from energy...")
    #energy = energy[energy.host < energy.neib]

    ### Removing entries where pair is not within cutoff distance ###
    LOGGER.info("remove non-neighbor cutoff pairs from energy...")
    energy = energy.merge(neibmap, on=["time", "host"])
    mask = []
    for res, nlist in zip(energy.neib, energy.nlist):
        mask.append(res in nlist)
    energy = energy[mask].drop(columns=["nlist"])
    del neibmap

    ### Frames okay? ###
    LOGGER.debug("neib:\n%s", neib)
    LOGGER.debug("order:\n%s", order)
    LOGGER.debug("energy:\n%s", energy)

    ### get dict to receive info for neighbor ###

    s_dict    = order.filter(["time", "host",  "host_Scd"])\
            .set_index(["time", "host"]).to_dict()["host_Scd"]
    resn_dict = order.filter(["time", "host", "host_type"])\
            .set_index(["time", "host"]).to_dict()["host_type"]

    ### Merging ###
    LOGGER.info("merging files...")

    order  = order.merge(neib, on=["time", "host", "host_type"])
    del neib
    final = energy.merge(order, on=["time", "host"])
    del order
    del energy
    final["neib_Scd"]  =  pd.Series(zip(final.time, final.neib)).map(s_dict)
    final["neib_type"] =  pd.Series(zip(final.time, final.neib)).map(resn_dict)

    final["DeltaScd"] = np.abs(final.host_Scd - final.neib_Scd)
    final["AvgScd"]   = ( final.host_Scd + final.neib_Scd ) / 2

    LOGGER.info("writing to csv...")

    for grpname, fr in final.groupby(["host_type", "neib_type"]):
        fr.to_csv( outputfile.format( '_'.join(grpname) ), index=False )


def create_nofs(outputfile="NofScd.csv",
    sfile="scd_distribution.dat",
    neighbortypefile="neighborcount.dat"):
    '''
        Merges order and neighbor data for each residue.
        Output columns are:
        <Time> <Host> <Lipid_type> <Host_Scd> <Ntot> <Ntype1> ... <NtypeX>

    '''
    def compare_cols(expectedcol, col):
        if set(col) != set(expectedcol):
            raise ValueError("Columns do not fit!\nexpected: {}\nfound:{}".format(expectedcol, col))

    LOGGER.info("loading order file (%s)...", sfile)
    #HEAD: time   resid   leaflet  resname   Scd
    order  = pd.read_table(sfile, sep="\s+")
    colnames = ["time", "resid", "leaflet", "resname", "Scd"]
    compare_cols(colnames, order.columns)

    LOGGER.info("loading neighbortype file (%s)...", neighbortypefile)
    #HEAD: time   resid resname   Neibs_Type1 Neibs_Type2 ...
    neib = pd.read_table(neighbortypefile, sep="\s+")
    if not ("time" in neib.columns and "resid" in neib.columns and "resname" in neib.columns):
        compare_cols(["time", "resid", "resname", "..."], neib.columns)

    ### Editing data ###
    neib["Ntot"] = neib.drop(columns=["time", "resid", "resname"]).sum(axis=1)

    ### Merging data ###
    final = neib.merge(order, on=["time", "resid", "resname"])

    ### Saving data ###
    final.to_csv(outputfile, index=False)
