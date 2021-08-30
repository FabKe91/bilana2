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
    energy = energy.rename(columns={"neighbor":"neib"})
    print(energy.info())

    LOGGER.info("loading neighbor info file (%s)...", neighborinfofile)
    #HEAD: resid  time  Number_of_neighbors List_of_Neighbors
    neibmap = pd.read_table(neighborinfofile, sep="\s+")
    colnames = ["time", "resid", "Number_of_neighbors", "List_of_Neighbors"]
    compare_cols(colnames, neibmap.columns)
    neibmap["nlist"] = neibmap.fillna('').List_of_Neighbors\
            .apply(lambda x: [int(i) for i in x.split(',') if i ])
    neibmap = neibmap.drop( columns=["Number_of_neighbors", "List_of_Neighbors"] )
    neibmap = neibmap.rename(columns={"resid":"host"})

    ### Removing entries where pair is not within cutoff distance ###
    energy = energy[energy.time >= neibmap.time.min()]
    energy_neibs = energy[["time", "host", "neib"]].merge(neibmap, on=["time", "host"])

    del neibmap

    LOGGER.info("remove non-neighbor cutoff pairs from energy...")
    mask = []
    for res, nlist in zip(energy_neibs.neib, energy_neibs.nlist):
        mask.append(res in nlist)
    
    del energy_neibs

    energy = energy[mask]
    print(energy.info())

    del mask


    ### Edit names | From here on work with new column list ! ###
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

    ### Frames okay? ###
    #LOGGER.debug("neib:\n%s", neib)
    #LOGGER.debug("order:\n%s", order)
    #LOGGER.debug("energy:\n%s", energy)

    ### get dict to receive info for neighbor ###
    LOGGER.info("loading order file (%s)...", sfile)
    #HEAD: time   resid   leaflet  resname   Scd
    order  = pd.read_table(sfile, sep="\s+")
    colnames = ["time", "resid", "leaflet", "resname", "Scd"]
    compare_cols(colnames, order.columns)
    order  = order.rename(columns={"resid":"host", "resname":"host_type", "Scd":"host_Scd"})

    s_dict    = order.filter(["time", "host",  "host_Scd"])\
            .set_index(["time", "host"]).to_dict()["host_Scd"]

    LOGGER.info("loading neighbortype file (%s)...", neighbortypefile)
    #HEAD: time   resid resname   Neibs_Type1 Neibs_Type2 ...
    neib = pd.read_table(neighbortypefile, sep="\s+")
    if not ("time" in neib.columns and "resid" in neib.columns and "resname" in neib.columns):
        compare_cols(["time", "resid", "resname", "..."], neib.columns)
    neib   = neib.rename(columns={"resid":"host", "resname":"host_type"})

    resn_dict = order.filter(["host", "host_type"])\
            .set_index("host").to_dict()["host_type"]

    ### Merging ###
    LOGGER.info("merging files...")

    order  = order.merge(neib, on=["time", "host", "host_type"])

    del neib

    final = energy.merge(order, on=["time", "host"])

    del energy

    final = final.merge(order.rename(columns={"host":"neib", "host_type":"neib_type", "host_Scd":"neib_Scd"}), on=["time", "neib"],suffixes=("_host", "_neib") ) # create neib_Scd and neib_type entries
    #final = final.merge(neib.rename(columns={"host":"neib", "host_type":"neib_type", }), on=["time", "neib", "neib_type"], ) # create neighborhood of neib entries 

    del order

    #final["neib_Scd"]  =  pd.Series(zip(final.time, final.neib)).map(s_dict)
    #final["neib_type"] =  pd.Series(zip(final.time, final.neib)).map(resn_dict)

    final["DeltaScd"] = np.abs(final.host_Scd - final.neib_Scd)
    final["AvgScd"]   = ( final.host_Scd + final.neib_Scd ) / 2

    ### get Ncb entry ###
    if "CHL1" in final.host_type.unique() or "CHL1" in final.neib_type.unique():
        final = final.merge(neibmap.rename(columns={"nlist":"host_nlist"}), on=["time", "host"])
        final = final.merge(neibmap.rename(columns={"host":"neib", "nlist":"neib_nlist"}), on=["time", "neib"])

        final["neibs"] = final.apply(lambda col: tuple(sorted(set(col["neib_nlist"] + col["host_nlist"]) - {col["host"]} - {col["neib"]} )), axis=1)
        final["neibs_both"] = final.apply(lambda col: tuple(sorted(set(col["neib_nlist"]) & set(col["host_nlist"]))), axis=1)

        final["neibl_type"]      = final.apply(lambda col: [resn_dict[i] for i in col["neibs"] ] , axis=1)
        final["neibl_type_both"] = final.apply(lambda col: [resn_dict[i] for i in col["neibs_both"] ] , axis=1)

        final["CHL1"] = final.apply(lambda col: col["neibl_type"].count("CHL1") , axis=1)
        final["CHL1_both"] = final.apply(lambda col: col["neibl_type_both"].count("CHL1") , axis=1)

        final = final.drop(columns=["host_nlist", "neib_nlist", "neibs", "neibs_both", "neibl_type", "neibl_type_both"])

    else:
        final["CHL1"] = 0
        final["CHL1_both"] = 0

    final["pair"] = np.array(['_'.join(i) for i in np.sort(final[["host_type", "neib_type"]].values, axis=1,)[:,::-1]])

    LOGGER.info("writing to csv...")
    print(final)

    for pairname, fr in final.groupby(["pair"]):
        fr.to_csv( outputfile.format( pairname ), index=False )


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
