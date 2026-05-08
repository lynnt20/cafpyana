from analysis_village.nuecc.makedf.make_nueccdf import * 

GRID_PARAMS = {
    "memory": "3GB",
    "cpu":    7,
    "disk":   "100GB",
    "lifetime": "1h",
}

DFS =   [make_nueccdf_base, make_hdrdf, make_potdf_bnb]
NAMES = ["nuecc", "hdr", "pot"]
