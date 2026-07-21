from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

# ============================================================================
# Helper functions
# ============================================================================

def NotInHigh(df): 
    fail = (df.z > 250) & (df.y > 100.) & (df.x < 0)
    return ~fail

def InFV_pfp(df):
    xmax = 195.
    zmin = 5.
    zmax = 490.
    ymax_highz = 100.
    pass_xz = (np.abs(df.x) < xmax) & (df.z > zmin) & (df.z < zmax)
    pass_y = ((df.z < 250) & (np.abs(df.y) < 195.)) | ((df.z > 250) & (df.y > -195.) & (df.y < ymax_highz))
    return pass_xz & pass_y

def slc_contained(pfpdf, slcdf, margin=0):
    # ignore the 'neutrino' pfp
    pfpdf = pfpdf[pfpdf.pfp.trk.producer!=4294967295]

    pfp_tpc0 = ((pfpdf.pfp.trk.start.x < -1*margin) & 
                (pfpdf.pfp.trk.end.x   < -1*margin) & 
                (pfpdf.pfp.shw.start.x < -1*margin) & 
                (pfpdf.pfp.shw.end.x   < -1*margin) & 
                (pfpdf.pfp.trk.start.x.notna()) &
                (pfpdf.pfp.trk.end.x.notna()) &
                (pfpdf.pfp.shw.start.x.notna()) &
                (pfpdf.pfp.shw.end.x.notna()))
    pfp_tpc1 = ((pfpdf.pfp.trk.start.x > margin) & 
                (pfpdf.pfp.trk.end.x   > margin) & 
                (pfpdf.pfp.shw.start.x > margin) & 
                (pfpdf.pfp.shw.end.x   > margin) & 
                (pfpdf.pfp.trk.start.x.notna()) &
                (pfpdf.pfp.trk.end.x.notna()) &
                (pfpdf.pfp.shw.start.x.notna()) &
                (pfpdf.pfp.shw.end.x.notna()))

    slc_tpc1 = pfp_tpc1.groupby(level=[0,1]).all()
    slc_tpc0 = pfp_tpc0.groupby(level=[0,1]).all()
    slc_cont = slc_tpc1 | slc_tpc0
    slcdf = multicol_add(slcdf, slc_cont.rename(("slc",'contained',f'margin_{margin}','tot')))
    slcdf = multicol_add(slcdf, slc_tpc0.rename(("slc",'contained',f'margin_{margin}','tpc0')))
    slcdf = multicol_add(slcdf, slc_tpc1.rename(("slc",'contained',f'margin_{margin}','tpc1')))
    return slcdf

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
    # keep only the maximum energy
    nu_df = nu_df.sort_values("E").groupby(level=[0,1]).nth(-1)
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

def make_mcnudf_nuecc_sig(f):
    mcdf = make_mcnudf_nuecc(f)
    signal_mask = ((InAV(mcdf.position, det="SBND")) &
                   (mcdf.iscc==1) &
                   (abs(mcdf.pdg)==12) &
                   (abs(mcdf.e.pdg)==11) )
    mcdf = mcdf[signal_mask]
    return mcdf

def make_mcnudf_nuecc_sigwgt(f, int_only=True,**kwargs):
    mcdf = make_mcnudf_nuecc(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    
    signal_mask = ((InFV(df=mcdf.position, inzback=0, det="SBND_nu26")) &
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
    DETECTOR = "SBND"
    # if (1 == det.unique()):
    #     DETECTOR = "SBND"
    # else:
    #     DETECTOR = "ICARUS"
    
    pfpdf = make_pfpdf(f)

    # load chi2pid for planes 0 and 1 for muon/proton (plane 2 already in trkbranches),
    # and all 3 planes for pion/kaon (not in trkbranches at all)
    chi2_extra = loadbranches(f["recTree"], [
        trkbranch + "chi2pid.0.chi2_muon",
        trkbranch + "chi2pid.0.chi2_proton",
        trkbranch + "chi2pid.0.chi2_pion",
        trkbranch + "chi2pid.0.chi2_kaon",
        trkbranch + "chi2pid.1.chi2_muon",
        trkbranch + "chi2pid.1.chi2_proton",
        trkbranch + "chi2pid.1.chi2_pion",
        trkbranch + "chi2pid.1.chi2_kaon",
        trkbranch + "chi2pid.2.chi2_pion",
        trkbranch + "chi2pid.2.chi2_kaon",
    ]).rec.slc.reco
    pfpdf = multicol_concat(pfpdf, chi2_extra)

    for _pid in ['muon', 'proton', 'pion', 'kaon']:
        _planes = [pfpdf[('pfp','trk','chi2pid',f'I{p}',f'chi2_{_pid}','')] for p in range(3)]
        _chi2_stack = pd.concat(_planes, axis=1)
        pfpdf[('pfp','trk','chi2pid','avg',f'chi2_{_pid}','')] = _chi2_stack.where(_chi2_stack > 0).mean(axis=1)

    slcdf = loadbranches(f["recTree"], slcbranches+barycenterFMbranches)
    slcdf = slcdf.rec

    slcdf = slc_contained(pfpdf, slcdf, margin=0)
    slcdf = slc_contained(pfpdf, slcdf, margin=5)
    slcdf = slc_contained(pfpdf, slcdf, margin=10)
    slcdf = slc_contained(pfpdf, slcdf, margin=20)
    slcdf = slc_contained(pfpdf, slcdf, margin=30)
    slcdf = slc_contained(pfpdf, slcdf, margin=50)
    slcdf = slc_contained(pfpdf, slcdf, margin=75)
    slcdf = slc_contained(pfpdf, slcdf, margin=100)
    
    pfp_NotInHigh = NotInHigh(pfpdf.pfp.trk.start) & NotInHigh(pfpdf.pfp.trk.end) & NotInHigh(pfpdf.pfp.shw.start) & NotInHigh(pfpdf.pfp.shw.end)
    slcdf = multicol_add(slcdf, pfp_NotInHigh.groupby(level=[0,1]).all().rename(('slc','pfp_notinhigh')))
    pfp_InFV = InFV_pfp(pfpdf.pfp.trk.start) & InFV_pfp(pfpdf.pfp.trk.end) & InFV_pfp(pfpdf.pfp.shw.start) & InFV_pfp(pfpdf.pfp.shw.end)
    slcdf = multicol_add(slcdf, pfp_InFV.groupby(level=[0,1]).all().rename(('slc','pfp_infv')))

    # get minimum abs x position of all pfps in the slice as a proxy for distance to the TPC boundary
    slcdf = multicol_add(slcdf, get_slcminx(pfpdf).rename(('slc','min_pfp_x')))

    pfpdf = pfpdf.drop('pfochar',axis=1,level=1)
    
    isshw = (pfpdf.pfp.trackScore < 0.5) & (pfpdf.pfp.shw.maxplane_energy > 0) & (pfpdf.pfp.trackScore > 0) & (pfpdf.pfp.shw.start.x == pfpdf.pfp.shw.start.x)
    istrk = (pfpdf.pfp.trackScore >= 0.5) & (pfpdf.pfp.trk.len > 0) & (pfpdf.pfp.trk.start.x == pfpdf.pfp.trk.start.x)
    isnon = ~(isshw | istrk)

    # min x of shower-type PFPs (signed, closest to TPC boundary)
    shw_sub = pfpdf[isshw]
    shw_xvals = np.array([shw_sub.pfp.shw.start.x.values, shw_sub.pfp.shw.end.x.values])
    shw_xvals = np.where(np.isnan(shw_xvals), np.inf, shw_xvals)
    shw_min_signed = np.choose(np.argmin(np.abs(shw_xvals), axis=0), shw_xvals)
    shw_min_x = pd.Series(shw_min_signed, index=shw_sub.index).groupby(level=[0,1]).apply(
        lambda g: g.iloc[np.abs(g).argmin()]
    )
    slcdf = multicol_add(slcdf, shw_min_x.rename(('slc', 'min_shw_x')))

    # min x of track-type PFPs (signed, closest to TPC boundary)
    trk_sub = pfpdf[istrk]
    trk_xvals = np.array([trk_sub.pfp.trk.start.x.values, trk_sub.pfp.trk.end.x.values])
    trk_xvals = np.where(np.isnan(trk_xvals), np.inf, trk_xvals)
    trk_min_signed = np.choose(np.argmin(np.abs(trk_xvals), axis=0), trk_xvals)
    trk_min_x = pd.Series(trk_min_signed, index=trk_sub.index).groupby(level=[0,1]).apply(
        lambda g: g.iloc[np.abs(g).argmin()]
    )
    slcdf = multicol_add(slcdf, trk_min_x.rename(('slc', 'min_trk_x')))

    # all classified (shw/trk) PFPs in the active volume
    shw_inav = (InAV(pfpdf[isshw].pfp.shw.start, det=DETECTOR) &
                InAV(pfpdf[isshw].pfp.shw.end,   det=DETECTOR))
    trk_inav = (InAV(pfpdf[istrk].pfp.trk.start, det=DETECTOR) &
                InAV(pfpdf[istrk].pfp.trk.end,   det=DETECTOR))
    pfp_inav = pd.concat([shw_inav, trk_inav])
    slcdf = multicol_add(slcdf, pfp_inav.groupby(level=[0,1]).all().rename(('slc', 'pfp_inav')))

    # reco particle ID counts using plane-averaged chi2pid
    chi2_mu_avg = pfpdf[('pfp','trk','chi2pid','avg','chi2_muon','')]
    chi2_p_avg  = pfpdf[('pfp','trk','chi2pid','avg','chi2_proton','')]
    is_reco_mu = istrk & (chi2_mu_avg > 0) & (chi2_mu_avg < 25) & (chi2_p_avg > 100)
    is_reco_p  = istrk & ~is_reco_mu & (chi2_p_avg > 0) & (chi2_p_avg < 90)
    slcdf = multicol_add(slcdf, is_reco_mu.groupby(level=[0,1]).sum().rename(('slc', 'n_reco_mu')))
    slcdf = multicol_add(slcdf, is_reco_p.groupby(level=[0,1]).sum().rename(('slc', 'n_reco_p')))

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

def make_nueccdf(f):
    slcdf = make_nueccdf_base(f)
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.5]
    slcdf = slcdf[InFV(df=slcdf.slc.vertex, det="SBND_nu26", inzback=0)]    
    return slcdf

def make_nueccdf_threshold(f):
    slcdf = make_nueccdf(f)
    slcdf = slcdf[slcdf.primshw.shw.reco_energy > 0.5]
    return slcdf

def make_intime_opflash_mc(f):
    opflashdf = make_opflashdf(f)
    return opflashdf[(opflashdf.time > -10) & (opflashdf.time < 10)]

def make_intime_opflash_data(f):
    opflashdf = make_opflashdf(f)
    framedf = make_framedf(f)[['frameApplyAtCaf']]
    df = pd.merge(opflashdf.reset_index(),
                  framedf.reset_index(),
                  left_on='entry',
                  right_on='entry',
                  how="left")
    df.set_index(opflashdf.index.names, verify_integrity=True, inplace=True)
    df = df[(df.time > -10) & (df.time < 10)]
    return df

# ============================================================================
# Data functions
# ============================================================================

def _prepare_nueccdf_data(slcdf, f):
    slcdf = slcdf.drop('tmatch', axis=1, level=1)  # slc level
    slcdf = slcdf.drop('truth',  axis=1, level=2)  # pfp level

    framedf = make_framedf(f)[['frameApplyAtCaf']]
    df = multicol_merge(slcdf.reset_index(),
                        framedf.reset_index(),
                        left_on=[('entry', '', '', '', '', '')],
                        right_on=[('entry', '', '', '', '', '')],
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_nueccdf_base_data(f):
    return _prepare_nueccdf_data(make_nueccdf_base(f), f)

def make_nueccdf_data(f):
    return _prepare_nueccdf_data(make_nueccdf(f), f)

def make_nueccdf_threshold_data(f):
    return _prepare_nueccdf_data(make_nueccdf_threshold(f), f)

def make_nuecc_df_data_sideband(f):
    slcdf = make_nueccdf_threshold(f)
    slcdf = slcdf[slcdf.primshw.shw.conversion_gap > 2]
    slcdf = slcdf[(slcdf.primshw.shw.bestplane_dEdx>3) &
                  (slcdf.primshw.shw.bestplane_dEdx<6)]
    return _prepare_nueccdf_data(slcdf, f)

def make_nuecc_df_data_sideband_debug(f):
    slcdf = make_nueccdf_threshold(f)
    slcdf = slcdf[slcdf.primshw.shw.conversion_gap > 2]
    slcdf = slcdf[(slcdf.primshw.shw.bestplane_dEdx>3)]
    return _prepare_nueccdf_data(slcdf, f)

def make_nuecc_df_data_signal(f): 
    slcdf = make_nueccdf_threshold(f)
    slcdf = slcdf[(slcdf.primtrk.trk.len < 200) |
                  (slcdf.primtrk.trk.len.isna())]
    slcdf = slcdf[(slcdf.primshw.shw.conversion_gap < 2) &
                  (slcdf.primshw.shw.conversion_gap > 0.001)]
    slcdf = slcdf[(slcdf.primshw.shw.bestplane_dEdx > 1.25) &
                  (slcdf.primshw.shw.bestplane_dEdx < 2.5)]
    slcdf = slcdf[(slcdf.primshw.shw.open_angle > 0.03) & 
                  (slcdf.primshw.shw.open_angle < 0.15)]
    slcdf = slcdf[(slcdf.primshw.shw.len > 10) &
                  (slcdf.primshw.shw.len < 200)]
    return _prepare_nueccdf_data(slcdf, f)
    

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

def _add_weights_to_nueccdf(df, f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie", "g4"], ar23p=False):
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
    
    if "fsi" in wgt_types:
        fsiwgtdf = make_fsi_weight_df(f)
        if fsiwgtdf is not None and not fsiwgtdf.empty:
            # Load genie_evtrec_idx: maps (entry, nu_subentry) → GenieEvtRecTree entry
            genie_idx_raw = loadbranches(f["recTree"], ["rec.mc.nu.genie_evtrec_idx"])
            while genie_idx_raw.columns.nlevels > 1:
                genie_idx_raw.columns = genie_idx_raw.columns.droplevel(0)
            genie_idx_ser = genie_idx_raw.iloc[:, 0]  # Series: (entry, nu_sub) → genie_entry

            # Build (entry, nu_subentry) pairs for selected neutrinos
            # nu_indices has values=nu_subentry, index level 0 = CAF entry
            entry_vals  = nu_indices.index.get_level_values(0)
            nu_sub_vals = nu_indices.values
            lookup_keys = pd.MultiIndex.from_arrays([entry_vals, nu_sub_vals])

            genie_entries = genie_idx_ser.reindex(lookup_keys).fillna(-1).astype(int).values

            ha2025  = fsiwgtdf[("fsi", "hA2025",  "")].reindex(genie_entries).fillna(1.0).values
            ha2025c = fsiwgtdf[("fsi", "hA2025c", "")].reindex(genie_entries).fillna(1.0).values

            fsi_aligned = pd.DataFrame({
                ("fsi", "hA2025",  ""): ha2025,
                ("fsi", "hA2025c", ""): ha2025c,
            }, index=nu_indices.index)
            fsi_aligned.columns = pd.MultiIndex.from_tuples(fsi_aligned.columns)
            wgt_dfs.append(fsi_aligned)

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

def make_nueccdf_mc_wgt_fsi(f, multisim_nuniv=100, slim=False, **kwargs):
    """Base selection with MC truth, systematic weights, and hA2025 FSI reweight."""
    df = make_nueccdf_mc(f)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim,
                                   wgt_types=["bnb", "genie", "g4", "fsi"], **kwargs)

def make_nueccdf_threshold_mc_wgt(f, multisim_nuniv=100, slim=False, **kwargs):
    df = make_nueccdf_threshold_mc(f)
    return _add_weights_to_nueccdf(df, f, multisim_nuniv=multisim_nuniv, slim=slim, **kwargs)