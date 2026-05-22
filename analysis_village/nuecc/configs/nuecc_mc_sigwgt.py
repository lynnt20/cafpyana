from analysis_village.nuecc.makedf.make_nueccdf import * 
GRID_PARAMS = {
    "memory": "6GB",
    "cpu":    3,
    "disk":   "25GB",
    "lifetime": "1h",
}

DFS =   [make_hdrdf,
         make_nueccdf_threshold_mc_wgt,
         make_mcnudf_nuecc_sigwgt,
         make_mcnulite_df_nuecc

         ]

NAMES = ["hdr", 
         "nuecc",
         "mcnuecc"
         "nulite"
         ]
