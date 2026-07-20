from analysis_village.nuecc.makedf.make_nueccdf import * 

GRID_PARAMS = {
    "memory": "3GB",
    "cpu":    7,
    "disk":   "100GB",
    "lifetime": "1h",
}

DFS =   [make_nuecc_df_data_signal, make_hdrdf, make_potdf_bnb, make_triggerdf]
NAMES = ["nuecc", "hdr", "pot","trigger"]