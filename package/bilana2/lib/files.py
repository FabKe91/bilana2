import numpy as np
import pandas as pd

def create_eofs(outputfile="EofScd{}.csv",
    efile="all_energies.dat",
    sfile="scd_distribution.dat",
    neighbortypefile="neighborcount.dat"):
    ''' Merges all_energies and scd_distribution file with output columns:
        var outputfile must contain a format control symbol ({}), so the interaction between
        specific lipid types is separated.

    <time> <host> <host_Scd> <neib> <neib_Scd> <molparts> <DeltaScd> <AvgScd> <leaflet> \
         <host_type> <neib_type> <resname1> <resname2> ... <resnameX>

    '''

    #HEAD: time   host   neighbor  molparts     VdW       Coul       Etot
    energy = pd.read_table(efile, delim_whitespace=True)

    #HEAD: time   resid   leaflet  resname   Scd
    order  = pd.read_table(sfile, delim_whitespace=True)

    #HEAD: time   resid   Neibs_Type1 Neibs_Type2 ...
    neib = pd.read_table(neighbortypefile, delim_whitespace=True)

    ### Edit names and remove duplicate entries (E_12 == E_21)
    order  = order.merge(neib, on=["time", "resid"])
    order  = order.rename(columns={"resid":"host", "resname":"host_type", "Scd":"host_Scd"})
    energy = energy.rename(columns={"Neighbor":"Neib"})
    energy = energy[energy.host < energy.neib]

    s_dict    = order.filter(["time", "host",  "host_Scd"])\
            .set_index(["time", "host"]).to_dict()["host_Scd"]
    resn_dict = order.filter(["time", "host", "host_type"])\
            .set_index(["time", "host"]).to_dict()["host_type"]

    final = energy.merge(order, on=["time", "host"])
    final["neib_Scd"]  =  pd.Series(zip(final.time, final.neib)).map(s_dict)
    final["neib_type"] =  pd.Series(zip(final.time, final.neib)).map(resn_dict)

    final["DeltaScd"] = np.abs(final.host_Scd - final.neib_Scd)
    final["AvgScd"]   = ( final.host_Scd + final.neib_Scd ) / 2

    for grpname, fr in final.groupby(["host_type", "neib_type"]):
        fr.to_csv( outputfile.format( '_'.join(grpname) ), index=False )
