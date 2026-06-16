from analysis_village.nuecc.makedf.make_nueccdf import *
GRID_PARAMS = {
    "memory": "6GB",
    "cpu":    7,
    "disk":   "25GB",
    "lifetime": "2h",
}

DFS =   [make_mcnudf_nuecc_sig,
         make_nueccdf_mc_wgt_fsi,
         make_hdrdf,
         make_mcnulite_df_nuecc]
NAMES = ["mcnuecc",
         "nuecc",
         "hdr",
         "nulite"]
