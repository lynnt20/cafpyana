from makedf.makedf import *
from pyanalib.pandas_helpers import *
from makedf.util import *

def make_sbndtriggerdf(f):
    triggerdf = loadbranches(f["recTree"], sbndtriggerbranches).rec.hdr.triggerinfo
    return triggerdf

def make_beamopflashdf(f):
    opflashdf = loadbranches(f["recTree"], opflashbranches).rec.opflashes
    opflashdf = opflashdf[(opflashdf.time > -10) & (opflashdf.time < 10)]
    return opflashdf

def make_softtrigdf(f):
    soft = loadbranches(f["recTree"], sbndsofttrigbranches).rec.soft_trig
    return soft

def make_alltrigdf(f):
    triggerdf = make_sbndtriggerdf(f)
    softtrigdf = make_softtrigdf(f)
    print(softtrigdf.columns)
    print(triggerdf.columns)
    df = multicol_merge(triggerdf,softtrigdf,left_index=True,right_index=True,how='left',validate='one_to_one')
    return df

def make_potspilldf_bnb(f):
    pot = loadbranches(f["recTree"], bnbpotspillbranches).rec.hdr.spillbnbinfo
    return pot

def make_timingdf(f):
    triggerdf = make_sbndtriggerdf(f)
    potspilldf = make_potspilldf_bnb(f)
    softtrigdf = make_softtrigdf(f)
    timingdf = loadbranches(f["recTree"], sbndtimingbranches).rec.sbnd_timings
    df = multicol_merge(triggerdf,potspilldf,left_index=True,right_index=True)
    df = multicol_merge(df,softtrigdf,left_index=True,right_index=True)
    df = multicol_merge(df,timingdf,left_index=True,right_index=True)
    return df

def make_nudf_data(f):
    slcdf = make_nudf(f)
    # drop truth cols for data
    slcdf = slcdf.drop('tmatch', axis=1,level=1) # slc level

    ## keep the only relevant column (for now)
    framedf = make_framedf(f)[['frameApplyAtCaf']]
    
    df = multicol_merge(slcdf.reset_index(), 
                        framedf.reset_index(),
                        left_on=[('entry', '', '',)],
                        right_on=[('entry', '', '',)], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_nudf_mc(f):
    slcdf = make_nudf(f)
    mcdf = make_mcnudf(f)
    mcdf.columns = pd.MultiIndex.from_tuples([tuple(["slc", "truth"] + list(c)) for c in mcdf.columns])
    df = multicol_merge(slcdf.reset_index(), 
                        mcdf.reset_index(),
                        left_on=[('entry', '', '', '',), 
                                ('slc', 'tmatch', 'idx', '',)], 
                        right_on=[('entry', '', '', '',), 
                                ('rec.mc.nu..index', '', '', '')], 
                        how="left")
    df = df.set_index(slcdf.index.names, verify_integrity=True)
    return df

def make_nudf(f):
    det = loadbranches(f["recTree"], ["rec.hdr.det"]).rec.hdr.det
    # if (1 == det.unique()):
    #     DETECTOR = "SBND"
    # else:
    #     DETECTOR = "ICARUS"
    # assert DETECTOR == "SBND"
    DETECTOR="SBND"
    
    slcdf = loadbranches(f["recTree"], slcbranches+barycenterFMbranches)
    slcdf = slcdf.rec    
    # pre-selection cuts
    slcdf = slcdf[slcdf.slc.is_clear_cosmic==0]
    slcdf = slcdf[InFV(df=slcdf.slc.vertex, inzback=0, det=DETECTOR)]

    return slcdf 
