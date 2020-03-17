import numpy as np
import pandas as pd

def create_eofs(outputfile="EofScd{}.csv", efile="all_energies.dat", sfile="scd_distribution.dat"):
    ''' Merges all_energies and scd_distribution file with output columns:
        var outputfile must contain a format control symbol ({}), so the interaction between
        specific lipid types is separated.

    <Time> <Host> <Host_Scd> <Neib> <Neib_Scd> <Molparts> <DeltaScd> <AvgScd> <leaflet> \
         <Host_Type> <Neib_Type> <Neibs_Type1> <Neibs_Type2> ... <Neibs_TypeX>

    '''

    #HEAD: Time   Host   Neighbor  Molparts     VdW       Coul       Etot
    energy = pd.read_table(efile, delim_whitespace=True)

    #HEAD: Time   Residue   leaflet  Type   Scd   Neibs_Type1 Neibs_Type2 ...
    order  = pd.read_table(sfile, delim_whitespace=True)

    ### Edit names and remove duplicate entries (E_12 == E_21)
    order  = order.rename(columns={"Residue":"Host", "Type":"Host_Type", "Scd":"Host_Scd"})
    energy = energy.rename(columns={"Neighbor":"Neib"})
    energy = energy[energy.Host < energy.Neib]

    s_dict    = order.filter(["Time", "Host",  "Host_Scd"])\
            .set_index(["Time", "Host"]).to_dict()["Host_Scd"]
    resn_dict = order.filter(["Time", "Host", "Host_Type"])\
            .set_index(["Time", "Host"]).to_dict()["Host_Type"]

    final = energy.merge(order, on=["Time", "Host"])
    final["Neib_Scd"]  =  pd.Series(zip(final.Time, final.Neib)).map(s_dict)
    final["Neib_Type"] =  pd.Series(zip(final.Time, final.Neib)).map(resn_dict)

    final["DeltaScd"] = np.abs(final.Host_Scd - final.Neib_Scd)
    final["AvgScd"]   = ( final.Host_Scd + final.Neib_Scd ) / 2

    for grpname, fr in final.groupby(["Host_Type", "Neib_Type"]):
        fr.to_csv( outputfile.format( '_'.join(grpname) ), index=False )
