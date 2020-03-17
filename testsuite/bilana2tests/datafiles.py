from pkg_resources import resource_filename

CHOL_ORDER    = resource_filename(__name__, "data/dppc_chol20_300/scd_distribution_original.dat")
CHOL_ENERGY   = resource_filename(__name__, "data/dppc_chol20_300/all_energies_original.dat")
CHOL_NEIGHBOR = resource_filename(__name__, "data/dppc_chol20_300/neighbor_info_original")
CHOL_INPUT    = resource_filename(__name__, "data/dppc_chol20_300/inputfile")
CHOL_RESINDEX = resource_filename(__name__, "data/dppc_chol20_300/resindex_all_original.ndx")
CHOL_TILT     = resource_filename(__name__, "data/dppc_chol20_300/tilt_original.csv")
print(CHOL_TILT)

CHOL_TPR      = resource_filename(__name__, "data/mdfiles_test/tpr/dppc_chol20_300.tpr")
CHOL_TRJ      = resource_filename(__name__, "data/mdfiles_test/md_trj/dppc_chol20_300_whole.xtc")
CHOL_GRO      = resource_filename(__name__, "data/mdfiles_test/initial_coords/dppc_chol20.gro")

CHIM_ORDER    = resource_filename(__name__, "data/dppc_chim20_330/scd_distribution_original.dat")
CHIM_ENERGY   = resource_filename(__name__, "data/dppc_chim20_330/all_energies_original.dat")
CHIM_NEIGHBOR = resource_filename(__name__, "data/dppc_chim20_330/neighbor_info_original")
CHIM_INPUT    = resource_filename(__name__, "data/dppc_chim20_330/inputfile")

CHIM_TPR      = resource_filename(__name__, "data/mdfiles_test/tpr/dppc_chim20_330.tpr")
CHIM_TRJ      = resource_filename(__name__, "data/mdfiles_test/md_trj/dppc_chim20_330_whole.xtc")
CHIM_GRO      = resource_filename(__name__, "data/mdfiles_test/initial_coords/dppc_chim20.gro")

del resource_filename
