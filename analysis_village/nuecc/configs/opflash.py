from analysis_village.nuecc.makedf.make_nueccdf import * 

GRID_PARAMS = {
    "memory": "3GB",
    "cpu":    7,
    "disk":   "100GB",
    "lifetime": "1h",
}

DFS =   [make_hdrdf, 
         make_intime_opflash_mc,
         make_mcnulite_df_nuecc]
NAMES = ["hdr", "opflash","nulite"]