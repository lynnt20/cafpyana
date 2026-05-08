from analysis_village.nuecc.makedf.make_nueccdf import * 

GRID_PARAMS = {
    "memory": "5GB",
    "cpu":    7,
    "disk":   "25GB",
    "lifetime": "1h",
}

DFS =   [make_nueccdf_mc_wgt, make_hdrdf]
NAMES = ["nuecc", "hdr"]
