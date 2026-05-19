from analysis_village.nuecc.makedf.make_nueccdf import * 

GRID_PARAMS = {
    "memory": "3GB",
    "cpu":    7,
    "disk":   "100GB",
    "lifetime": "1h",
}

DFS =   [make_mcnudf_nuecc_sig, make_nueccdf_base_mc, make_hdrdf, make_potdf_bnb]
NAMES = ["mcnuecc","nuecc", "hdr", "pot"]
