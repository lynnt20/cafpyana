from analysis_village.nuecc.makedf.make_nueccdf import * 

DFS =   [make_nueccdf_withcuts_control_mc_wgt, 
         make_nueccdf_withcuts_mc_wgt,
         make_hdrdf, make_potdf_bnb, make_mcnudf_nuecc]
NAMES = ["nuecc_control", "nuecc",
         "hdr", "pot", "mcnu"]
