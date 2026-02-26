from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

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
    return mcdf


def make_sigmcnudf_wgt(f):
    mcdf = make_mcdf(f)
    mcdf["ind"] = mcdf.index.get_level_values(1)
    
    if 'mu'  in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('mu', axis=1,level=0)
    if 'p'   in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('p',  axis=1,level=0)
    if 'cpi' in list(zip(*list(mcdf.columns)))[0]:  mcdf = mcdf.drop('cpi',axis=1,level=0)
    
    ## select out signal events 
    mcdf = mcdf[InAV(df=mcdf.position)] # in AV
    mcdf = mcdf[(mcdf.iscc==1) & 
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

def make_nueccdf_mc_wgt(f):
    df = make_nueccdf_mc(f,include_weights=True)
    return df

def make_nueccdf_mc_wgt_postsel(f, multisim_nuniv=100, slim=False, wgt_types=["bnb", "genie"], **kwargs):
    """
    Generates weights only for neutrino interactions that pass preselection.
    """
    # First, build the mc df WITHOUT weights
    df = make_nueccdf_mc(f, include_weights=False)
    
    # Get only the unique truth-matched neutrino indices that survived preselection
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

def make_nueccdf_mc(f, include_weights=False,multisim_nuniv=100,slim=False,**kwargs):
    
    slcdf = make_nueccdf(f)
    mcdf = make_mcnudf_nuecc(f,include_weights=include_weights,multisim_nuniv=multisim_nuniv,slim=slim,**kwargs)
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

def make_nueccdf_mc_withcuts(f):
    slcdf = make_nueccdf(f)
    slcdf = slcdf[np.isnan(slcdf.primtrk.trk.len) | (slcdf.primtrk.trk.len < 200)]
    slcdf = slcdf.drop('primtrk',axis=1,level=0)
    slcdf = slcdf[slcdf.primshw.shw.conversion_gap < 2]
    slcdf = slcdf[(slcdf.primshw.shw.bestplane_dEdx > 1) & (slcdf.primshw.shw.bestplane_dEdx < 2.5)]
    slcdf = slcdf[slcdf.primshw.shw.open_angle < 0.2] 
    
    mcdf = make_mcnudf_nuecc(f,include_weights=False)
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
    
    # pre-selection cuts
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[slcdf.slc.nu_score > 0.4]
    slcdf = slcdf[InAV(df=slcdf.slc.vertex, det=DETECTOR)]    
    
    return slcdf
