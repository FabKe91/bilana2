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
    order  = order.rename(columns={"resid":"Host", "resname":"Host_Type", "Scd":"Host_Scd"})
    energy = energy.rename(columns={"Neighbor":"Neib"})
    energy = energy[energy.Host < energy.Neib]

    s_dict    = order.filter(["time", "host",  "host_Scd"])\
            .set_index(["time", "host"]).to_dict()["host_Scd"]
    resn_dict = order.filter(["time", "host", "host_Type"])\
            .set_index(["time", "host"]).to_dict()["host_Type"]

    final = energy.merge(order, on=["time", "host"])
    final["neib_Scd"]  =  pd.Series(zip(final.Time, final.Neib)).map(s_dict)
    final["neib_type"] =  pd.Series(zip(final.Time, final.Neib)).map(resn_dict)

    final["DeltaScd"] = np.abs(final.Host_Scd - final.Neib_Scd)
    final["AvgScd"]   = ( final.Host_Scd + final.Neib_Scd ) / 2

    for grpname, fr in final.groupby(["host_type", "neib_type"]):
        fr.to_csv( outputfile.format( '_'.join(grpname) ), index=False )
