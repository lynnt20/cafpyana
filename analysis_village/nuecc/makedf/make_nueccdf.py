from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

# ============================================================================
# Helper functions
# ============================================================================

def get_pfpcontained(pfpdf, margin=0.0):
    """
    Compute pfpcontained flag for slices with optional distance margin.

    Checks if all PFPs in a slice satisfy containment criteria. A PFP is
    considered contained if both its track and shower don't cross the
    specified boundary.

    Parameters
    ----------
    pfpdf : pandas.DataFrame
        PFP-level dataframe with MultiIndex [run, event, pfp]
    margin : float, optional
        Distance margin in cm from YZ=0 plane (default: 0.0).
        E.g., margin=5 requires |x| > 5 for both start and end points.

    Returns
    -------
    pandas.Series
        Boolean series indexed by [run, event] indicating if all pfps
        in the slice are contained
    """

    trk_contained = (pfpdf.pfp.trk.start.x > 0) == (pfpdf.pfp.trk.end.x > 0)
    shw_contained = (pfpdf.pfp.shw.start.x > 0) == (pfpdf.pfp.shw.end.x > 0)
    if margin == 0.0:
        return (trk_contained & shw_contained).groupby(level=[0, 1]).all()
    else:
        trk_margin = (abs(pfpdf.pfp.trk.start.x) > margin) == (abs(pfpdf.pfp.trk.end.x) > margin)
        shw_margin = (abs(pfpdf.pfp.shw.start.x) > margin) == (abs(pfpdf.pfp.shw.end.x) > margin)
        return (trk_contained & shw_contained & trk_margin & shw_margin).groupby(level=[0, 1]).all()

def get_slcminx(pfpdf):
    # get minimum x position of all pfps in the slice as a proxy for distance to the TPC boundary
    # but keep the sign to distinguish proximity to upstream vs downstream boundary
    endpoints = np.array([
        pfpdf.pfp.trk.start.x,
        pfpdf.pfp.trk.end.x,
        pfpdf.pfp.shw.start.x,
        pfpdf.pfp.shw.end.x,
    ])
    endpoints = np.where(np.isnan(endpoints), np.inf, endpoints)  # treat NaNs as infinitely far from boundary
    endpoints_abs = np.abs(endpoints)
    min_abs_per_pfp = np.min(endpoints_abs, axis=0)
    min_idx_per_pfp = np.argmin(endpoints_abs, axis=0)
    min_signed_per_pfp = np.choose(min_idx_per_pfp, endpoints)
    
    min_signed_series = pd.Series(min_signed_per_pfp, index=pfpdf.index)
    result = min_signed_series.groupby(level=[0,1]).apply(
        lambda group: group.iloc[np.abs(group).argmin()]
    )
    return result

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
    # if 'mu'  in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('mu', axis=1,level=0)
    if 'p'   in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('p',  axis=1,level=0)
    if 'cpi' in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('cpi',axis=1,level=0)
    
    mcdf.loc[:, ('e','totp','')] = np.sqrt(mcdf.e.genp.x**2 + mcdf.e.genp.y**2 + mcdf.e.genp.z**2)
    mcdf.loc[:, ('e','dir','x')] = mcdf.e.genp.x/mcdf.e.totp
    mcdf.loc[:, ('e','dir','y')] = mcdf.e.genp.y/mcdf.e.totp
    mcdf.loc[:, ('e','dir','z')] = mcdf.e.genp.z/mcdf.e.totp
    
    mcdf.loc[:, ('pi0','totp','')] = np.sqrt(mcdf.pi0.genp.x**2 + mcdf.pi0.genp.y**2 + mcdf.pi0.genp.z**2)
    mcdf.loc[:, ('pi0','dir','x')] = mcdf.pi0.genp.x/mcdf.pi0.totp
    mcdf.loc[:, ('pi0','dir','y')] = mcdf.pi0.genp.y/mcdf.pi0.totp
    mcdf.loc[:, ('pi0','dir','z')] = mcdf.pi0.genp.z/mcdf.pi0.totp
    
    return mcdf

def make_mcnudf_nuecc_sigwgt(f, int_only=True,**kwargs):
    mcdf = make_mcnudf_nuecc(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    
    signal_mask = ((InFV(df=mcdf.position, inzback=0, det="SBND_nohighyz")) &
                   (mcdf.iscc==1) &
                   (abs(mcdf.pdg)==12) &
                   (abs(mcdf.e.pdg)==11) &
                   (mcdf.e.genE > 0.5)) # nueCC signal definition
    if int_only:
        mcdf = mcdf[signal_mask]
    else:
        # select first-level groups that contain at least one signal event
        sig_groups = mcdf[signal_mask].index.get_level_values(0).unique()
        mcdf = mcdf[mcdf.index.get_level_values(0).isin(sig_groups)]
    
    if mcdf.empty:
        return mcdf  # skip weight merging if no signal events
    geniewgtdf = geniesyst.geniesyst(f, 
                                     mcdf.ind, 
                                     multisim_nuniv=100, 
                                     slim=False, 
                                     systematics=None,
                                     **kwargs)
    mcdf = multicol_concat(mcdf, geniewgtdf)
    return mcdf

def make_mcnudf_nuecc_sigwgt_ar23p(f):
    # get ar23p weights for all interactions with a signal event 
    return make_mcnudf_nuecc_sigwgt(f, int_only=True,ar23p=True)

def make_mcnudf_nuecc_sigwgt_ar23p_only(f):
    # get ar23p weights for all interactions with a signal event 
    return make_mcnudf_nuecc_sigwgt(f, int_only=False,ar23p_only=True)

# ============================================================================
# Base selection functions (call hierarchy)
# ============================================================================

def make_nueccdf_base(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    if (1 == det.unique()):
        DETECTOR = "SBND"
    else:
        DETECTOR = "ICARUS"
    
    pfpdf = make_pfpdf(f)

    slcdf = loadbranches(f["recTree"], slcbranches+barycenterFMbranches)
    slcdf = slcdf.rec

    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=0.0).rename(('slc','contained','0cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=5.0).rename(('slc','contained','5cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=10.0).rename(('slc','contained','10cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=20.0).rename(('slc','contained','20cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=30.0).rename(('slc','contained','30cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=50.0).rename(('slc','contained','50cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=75.0).rename(('slc','contained','75cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=100.0).rename(('slc','contained','100cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=125.0).rename(('slc','contained','125cm')))
    slcdf = multicol_add(slcdf, get_pfpcontained(pfpdf, margin=150.0).rename(('slc','contained','150cm')))
    
    # get minimum abs x position of all pfps in the slice as a proxy for distance to the TPC boundary
    slcdf = multicol_add(slcdf, get_slcminx(pfpdf).rename(('slc','min_pfp_x')))

    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    
    isshw = (pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0) & (pfpdf.pfp.trackScore > 0) & (pfpdf.pfp.shw.start.x == pfpdf.pfp.shw.start.x) 
    istrk = (pfpdf.pfp.trackScore >= 0.5) & (pfpdf.pfp.trk.len > 0) & (pfpdf.pfp.trk.start.x == pfpdf.pfp.trk.start.x)
    isnon = ~(isshw | istrk) 

    slcdf = multicol_add(slcdf, isshw.groupby(level=[0,1]).sum().rename(('slc','nshw')))
    slcdf = multicol_add(slcdf, istrk.groupby(level=[0,1]).sum().rename(('slc','ntrk')))
    slcdf = multicol_add(slcdf, isnon.groupby(level=[0,1]).sum().rename(('slc','nnon')))
    
    slcdf = multicol_add(slcdf, pfpdf[isshw].pfp.shw.maxplane_energy.groupby(level=[0,1]).sum().rename(('slc','tot_shw_energy')))
    
    ## primary shw candidate is shw pfp with highest energy, valid energy, and score < 0.5
    shwdf = pfpdf[isshw].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from trk attributes
    shwdf = shwdf.drop('trk',axis=1,level=1)
    shwdf.columns = shwdf.columns.set_levels(['primshw'],level=0)
    slcdf = multicol_merge(slcdf, shwdf.droplevel(-1),left_index=True,right_index=True,how="right",validate="one_to_one")

    ## secondary shower is shw pfp with second highest energy, valid energy, and score < 0.5 
    shwsecdf = pfpdf[isshw].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','shw','maxplane_energy','','','')]).groupby(level=[0,1]).nth(-2)
    shwsecdf = shwsecdf.drop('trk',axis=1,level=1)
    shwsecdf.columns = shwsecdf.columns.set_levels(['secshw'],level=0)
    slcdf = multicol_merge(slcdf, shwsecdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    ## primary trk is track pfp with the longest length
    trkdf = pfpdf[istrk].sort_values(pfpdf.pfp.index.names[:-1] + [('pfp','trk','len','','','')]).groupby(level=[0,1]).nth(-1)
    # drop all columns that are from shw attributes
    trkdf = trkdf.drop('shw',axis=1,level=1)
    trkdf.columns = trkdf.columns.set_levels(['primtrk'],level=0)
    slcdf = multicol_merge(slcdf, trkdf.droplevel(-1),left_index=True,right_index=True,how="left",validate="one_to_one")

    # add a shower energy variable that applies a scale factor to the max plane energy of the primary shower candidate
    shower_scale=1.17
    slcdf = multicol_add(slcdf,(slcdf.primshw.shw.maxplane_energy*shower_scale).rename(("primshw","shw","reco_energy")))
    
    return slcdf

# ============================================================================
# Data functions
# ============================================================================

def make_nueccdf(f):
    slcdf = make_nueccdf_base(f)
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]
    slcdf = slcdf[InFV(df=slcdf.slc.vertex, det="SBND_nohighyz", inzback=0)]    
    return slcdf

def make_nueccdf_threshold(f):
    slcdf = make_nueccdf(f)
    slcdf = slcdf[slcdf.primshw.shw.reco_energy > 0.5]
    return slcdf
    
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

def make_nueccdf_base_mc(f):
    slcdf = make_nueccdf_base(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f)

def make_nueccdf_mc(f):
    slcdf = make_nueccdf(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f)

def make_nueccdf_threshold_mc(f):
    slcdf = make_nueccdf_threshold(f)
    return _merge_nueccdf_with_mc_truth(slcdf, f)

# ============================================================================
# Systematic weights helper
# ============================================================================

def _add_weights_to_nueccdf(df, f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie","g4"],ar23p=False):
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

    # No selected neutrinos in this file: skip weight extraction entirely.
    if nu_indices.empty:
        return df
    
    wgt_dfs = []
    
    if "genie" in wgt_types:
        geniewgtdf = geniesyst.geniesyst(f, 
                                         nu_indices, 
                                         multisim_nuniv=multisim_nuniv, 
                                         slim=slim, 
                                         systematics=None,
                                         ar23p=ar23p)
        if geniewgtdf is not None and geniewgtdf.shape[1] > 0:
            wgt_dfs.append(geniewgtdf)
    
    if "bnb" in wgt_types:
        bnbwgtdf = bnbsyst.bnbsyst(f, 
                                    nu_indices, 
                                    multisim_nuniv=multisim_nuniv, 
                                    slim=slim)
        if bnbwgtdf is not None and bnbwgtdf.shape[1] > 0:
            wgt_dfs.append(bnbwgtdf)
    
    if "g4" in wgt_types:
        g4wgtdf = g4syst.g4syst(f, 
                                 nu_indices, 
                                 multisim_nuniv=multisim_nuniv, 
                                 slim=slim)
        if g4wgtdf is not None and g4wgtdf.shape[1] > 0:
            wgt_dfs.append(g4wgtdf)
    
    if wgt_dfs:
        wgtdf = pd.concat(wgt_dfs, axis=1)
        del wgt_dfs  # free intermediate list
        if wgtdf.shape[1] == 0:
            return df
        wgtdf.columns = pd.MultiIndex.from_tuples(
            [tuple(["slc", "truth"] + list(c)) for c in wgtdf.columns]
        )
        df = multicol_concat(df, wgtdf)
        del wgtdf
    
    return df

# ============================================================================
# MC functions with weights
# ============================================================================

def make_nueccdf_mc_wgt(f, multisim_nuniv=100, slim=False, **kwargs):
    """
    Base selection with MC truth and systematic weights.
    Weights are calculated for selected indices only to reduce overhead.
    """
    df = make_nueccdf_mc(f)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim, **kwargs)

def make_nueccdf_mc_wgt_ar23(f, multisim_nuniv=100, slim=False, **kwargs):
    """
    Base selection with MC truth and systematic weights.
    Weights are calculated for selected indices only to reduce overhead.
    """
    df = make_nueccdf_mc(f)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim, ar23p=True,**kwargs)