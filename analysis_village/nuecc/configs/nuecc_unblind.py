from analysis_village.nuecc.makedf.make_nueccdf import * 

GRID_PARAMS = {
    "memory": "3GB",
    "cpu":    7,
    "disk":   "100GB",
    "lifetime": "1h",
}

DFS =   [make_nueccdf_threshold_data, make_hdrdf, make_potdf_bnb, make_triggerdf, make_intime_opflash_data]
NAMES = ["nuecc", "hdr", "pot","trigger", "opflash"]