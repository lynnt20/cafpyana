from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

# ============================================================================
# Independent base functions
# ============================================================================

def make_mcnulite_df_nuecc(f):
    rse_df = make_hdrdf(f)
    nu_df = loadbranches(f["recTree"], ["rec.mc.nu.E","rec.mc.nu.pdg"]).rec.mc.nu
    names = nu_df.index.names
    df = nu_df.reset_index().merge(rse_df.reset_index(),on='entry').set_index(names)
    return df 

def make_mcnudf_nuecc(f,**args):
    mcdf = make_mcnudf(f,**args)
    # drop mcdf columns not relevant for this analysis
    if 'mu'  in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('mu', axis=1,level=0)
    if 'p'   in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('p',  axis=1,level=0)
    if 'cpi' in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('cpi',axis=1,level=0)
    
    mcdf.loc[:, ('e','totp','')] = np.sqrt(mcdf.e.genp.x**2 + mcdf.e.genp.y**2 + mcdf.e.genp.z**2)

    # opening angles
    mcdf.loc[:, ('e','dir','x')] = mcdf.e.genp.x/mcdf.e.totp
    mcdf.loc[:, ('e','dir','y')] = mcdf.e.genp.y/mcdf.e.totp
    mcdf.loc[:, ('e','dir','z')] = mcdf.e.genp.z/mcdf.e.totp
    return mcdf

def make_mcnudf_nuecc_sigwgt(f):
    mcdf = make_mcnudf_nuecc(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    ## select out signal events 
    mcdf = mcdf[(InFV(df=mcdf.position, inzback=0, det="SBND_nohighyz")) &
                (mcdf.iscc==1) & 
                (abs(mcdf.pdg)==12) & 
                (abs(mcdf.e.pdg)==11) & 
                (mcdf.e.genE > 0.5)] # nueCC signal definition
    
    geniewgtdf = geniesyst.geniesyst(f, 
                                     mcdf.ind, 
                                     multisim_nuniv=100, 
                                     slim=False, 
                                     systematics=None)
    mcdf = multicol_concat(mcdf, geniewgtdf)
    return mcdf

# ============================================================================
# Base selection functions (call hierarchy)
# ============================================================================

def make_nueccdf_debug(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    pfpdf = make_pfpdf(f)
    slcdf = loadbranches(f["recTree"], slcbranches+barycenterFMbranches)
    slcdf = slcdf.rec
    
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    ## primary shw candidate is shw pfp with highest energy, valid energy, and score < 0.5
    shwdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from trk attributes
    shwdf = shwdf.drop('trk',axis=1,level=1)
    shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
    slcdf = multicol_merge(slcdf, shwdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    ## primary trk is track pfp with the longest length
    trkdf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','trk','len','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    trkdf = trkdf.drop('shw',axis=1,level=1)
    trkdf.columns = trkdf.columns.set_levels(['primtrk'],level=0)
    slcdf = multicol_merge(slcdf, trkdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    ## secondary shower is shw pfp with second highest energy, valid energy, and score < 0.5 
    shwsecdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-2)
    shwsecdf = shwsecdf.drop('trk',axis=1,level=1)
    shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'],level=0)
    slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
    
    return slcdf

def make_nueccdf(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"

    assert DETECTOR == "SBND"
    
    pfpdf = make_pfpdf(f)
    slcdf = loadbranches(f["recTree"], slcbranches+barycenterFMbranches)
    slcdf = slcdf.rec
    
    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    ## primary shw candidate is shw pfp with highest energy, valid energy, and score < 0.5
    shwdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from trk attributes
    shwdf = shwdf.drop('trk',axis=1,level=1)
    shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
    slcdf = multicol_merge(slcdf, shwdf.droplevel(-1),left_index=True,right_index=True,how="right",validate="one_to_one")

    ## primary trk is track pfp with the longest length
    trkdf = pfpdf[(pfpdf.pfp.trackScore > 0.5) & (pfpdf.pfp.trk.len > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','trk','len','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    trkdf = trkdf.drop('shw',axis=1,level=1)
    trkdf.columns = trkdf.columns.set_levels(['primtrk'],level=0)
    slcdf = multicol_merge(slcdf, trkdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    ## secondary shower is shw pfp with second highest energy, valid energy, and score < 0.5 
    shwsecdf = pfpdf[(pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0)].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-2)
    shwsecdf = shwsecdf.drop('trk',axis=1,level=1)
    shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'],level=0)
    slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")
    
    # add a shower energy variable that applies a scale factor to the max plane energy of the primary shower candidate
    shower_scale=1.25
    slcdf = multicol_add(slcdf,(slcdf.primshw.shw.maxplane_energy*shower_scale).rename(("primshw","shw","reco_energy")))
    
    # pre-selection cuts
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.4]
    slcdf = slcdf[InFV(df=slcdf.slc.vertex, det="SBND_nohighyz", inzback=0)]    
    
    return slcdf

def make_nueccdf_withcuts(f):

    score_cut=0.02
    min_shower_energy=0.5
    
    max_track_length=200
    min_conversion_gap=0.001
    max_conversion_gap=2
    min_dedx=1.25
    max_dedx=2.5
    
    df = make_nueccdf(f)
    df = df[df.slc.barycenterFM.score > score_cut]
    
    df = df[np.isnan(df.primtrk.trk.len) | (df.primtrk.trk.len < max_track_length)]

    df = df[(df.primshw.shw.conversion_gap < max_conversion_gap) & 
            (df.primshw.shw.conversion_gap > min_conversion_gap)]

    df = df[(df.primshw.shw.bestplane_dEdx > min_dedx) & (df.primshw.shw.bestplane_dEdx < max_dedx)]
    
    return df

def make_nueccdf_withcuts_control(f):

    score_cut=0.02
    min_shower_energy=0.5
    
    max_track_length=1e5
    min_conversion_gap=2
    max_conversion_gap=1e9
    min_dedx=3
    max_dedx=6
    
    df = make_nueccdf(f)
    df = df[df.slc.barycenterFM.score > score_cut]
    
    df = df[np.isnan(df.primtrk.trk.len) | (df.primtrk.trk.len < max_track_length)]

    df = df[(df.primshw.shw.conversion_gap < max_conversion_gap) & 
            (df.primshw.shw.conversion_gap > min_conversion_gap)]

    df = df[(df.primshw.shw.bestplane_dEdx > min_dedx) & (df.primshw.shw.bestplane_dEdx < max_dedx)]
    
    return df

# ============================================================================
# Data functions
# ============================================================================

def make_nueccdf_data(f):
    slcdf = make_nueccdf(f)
    # drop truth cols for data
    slcdf = slcdf.drop('tmatch', axis=1,level=1) # slc level
    slcdf = slcdf.drop('truth',  axis=1,level=2) # pfp level
    
    ## keep the only relevant column (for now)
    framedf = make_framedf(f)[['frameApplyAtCaf']]
    
    df = multicol_merge(slcdf.reset_index(), 
                        framedf.reset_index(),
                        left_on=[('entry', '', '', '', '', '')],
                        right_on=[('entry', '', '', '', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_nueccdf_debug_data(f):
    slcdf = make_nueccdf_debug(f)
    # drop truth cols for data
    slcdf = slcdf.drop('tmatch', axis=1,level=1) # slc level
    slcdf = slcdf.drop('truth',  axis=1,level=2) # pfp level
    
    ## keep the only relevant column (for now)
    framedf = make_framedf(f)[['frameApplyAtCaf']]
    
    df = multicol_merge(slcdf.reset_index(), 
                        framedf.reset_index(),
                        left_on=[('entry', '', '', '', '', '')],
                        right_on=[('entry', '', '', '', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

# ============================================================================
# MC truth merge helper
# ============================================================================

def _merge_nueccdf_with_mc_truth(slcdf, f, include_weights=False, multisim_nuniv=100, slim=False, **kwargs):
    """
    Helper function to merge a slc dataframe with MC truth information.
    
    Parameters
    ----------
    slcdf : pandas.DataFrame
        Slice-level dataframe to merge truth information into
    f : ROOT file
        Input ROOT file
    include_weights : bool, optional
        Whether to include weights (default: False)
    multisim_nuniv : int, optional
        Number of multisim universes (default: 100)
    slim : bool, optional
        Whether to slim the output (default: False)
    **kwargs : dict
        Additional keyword arguments passed to make_mcnudf_nuecc
    
    Returns
    -------
    pandas.DataFrame
        Merged dataframe with truth information
    """
    mcdf = make_mcnudf_nuecc(f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, slim=slim, **kwargs)
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "truth"] + list(c)) for c in mcdf.columns])
    df = multicol_merge(slcdf.reset_index(), 
                        mcdf.reset_index(),
                        left_on=[('entry', '', '', '', '', ''), 
                                ('slc', 'tmatch', 'idx', '', '', '')], 
                        right_on=[('entry', '', '', '', '', ''), 
                                ('rec.mc.nu..index', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

# ============================================================================
# MC functions (with truth matching)
# ============================================================================

def make_nueccdf_mc(f, include_weights=False, multisim_nuniv=100, slim=False, **kwargs):
    """
    Merge base selection with MC truth information.
    """
    slcdf = make_nueccdf(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, slim=slim, **kwargs)

def make_nueccdf_debug_mc(f, include_weights=False, multisim_nuniv=100, slim=False, **kwargs):
    """
    Merge base selection with MC truth information for debug version of df.
    """
    slcdf = make_nueccdf_debug(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f, include_weights=include_weights, multisim_nuniv=multisim_nuniv, slim=slim, **kwargs)

def make_nueccdf_withcuts_mc(f):
    """
    Merge pre-selected signal region with MC truth information.
    """
    slcdf = make_nueccdf_withcuts(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f, include_weights=False)

def make_nueccdf_withcuts_control_mc(f):
    """
    Merge pre-selected control region with MC truth information.
    """
    slcdf = make_nueccdf_withcuts_control(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f, include_weights=False)

# ============================================================================
# Systematic weights helper
# ============================================================================

def _add_weights_to_nueccdf(df, f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie"]):
    """
    Helper function to add systematic weights to a neutrino CC DataFrame.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe with truth-matched information
    f : ROOT file
        Input ROOT file
    multisim_nuniv : int, optional
        Number of multisim universes (default: 100)
    slim : bool, optional
        Whether to slim the output (default: False)
    wgt_types : list, optional
        List of weight types to include, e.g. ["bnb", "genie"] (default: ["bnb", "genie"])
    
    Returns
    -------
    pandas.DataFrame
        DataFrame with added weight columns
    """
    # Get only the unique truth-matched neutrino indices that survived selection
    nu_indices = df[('slc', 'truth', 'ind', '', '', '')].dropna().astype(int)
    nu_indices = nu_indices[~nu_indices.index.duplicated()]  # deduplicate
    
    wgt_dfs = []
    
    if "genie" in wgt_types:
        geniewgtdf = geniesyst.geniesyst(f, 
                                         nu_indices, 
                                         multisim_nuniv=multisim_nuniv, 
                                         slim=slim, 
                                         systematics=None)
        wgt_dfs.append(geniewgtdf)
    
    if "bnb" in wgt_types:
        bnbwgtdf = bnbsyst.bnbsyst(f, 
                                    nu_indices, 
                                    multisim_nuniv=multisim_nuniv, 
                                    slim=slim)
        wgt_dfs.append(bnbwgtdf)
    
    if wgt_dfs:
        wgtdf = pd.concat(wgt_dfs, axis=1)
        del wgt_dfs  # free intermediate list
        wgtdf.columns = pd.MultiIndex.from_tuples(
            [tuple(["slc", "truth"] + list(c)) for c in wgtdf.columns]
        )
        df = multicol_concat(df, wgtdf)
        del wgtdf
    
    return df

# ============================================================================
# MC functions with weights
# ============================================================================

def make_nueccdf_mc_wgt(f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie"]):
    """
    Base selection with MC truth and systematic weights.
    Weights are calculated for selected indices only to reduce overhead.
    """
    df = make_nueccdf_mc(f, include_weights=False)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim, wgt_types=wgt_types)

def make_nueccdf_withcuts_mc_wgt(f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie"], **kwargs):
    """
    Pre-selected signal region with MC truth and systematic weights.
    Weights are calculated for selected indices only to reduce overhead.
    """
    df = make_nueccdf_withcuts_mc(f)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim, wgt_types=wgt_types)

def make_nueccdf_withcuts_control_mc_wgt(f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie"], **kwargs):
    """
    Pre-selected control region with MC truth and systematic weights.
    Weights are calculated for selected indices only to reduce overhead.
    """
    df = make_nueccdf_withcuts_control_mc(f)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim, wgt_types=wgt_types)