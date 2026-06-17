#!/usr/bin/env python3
"""
Compute volume-averaged BNB gsimple flux for defined FV volumes and save to pickle.

Output pickle contains:
  spectra    : dict[vname][flavor] -> np.ndarray, shape (N_E_BINS,)
               flux in m^-2 POT^-1 per 50 MeV bin
  bin_edges  : np.ndarray, shape (N_E_BINS+1,)
  bin_centers: np.ndarray, shape (N_E_BINS,)
  total_pot  : float
  flavors    : list[str]
  volumes    : dict[vname] -> list of box dicts
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing
import os
import pickle

import numpy as np
import uproot
from tqdm import tqdm

from raytrace_volume_defs import RAYTRACE_VOLUME_DEFS


def path_length_box(vx, vy, vz, dx, dy, dz, box):
    def slab_t(v, d, lo, hi):
        with np.errstate(divide="ignore", invalid="ignore"):
            t1 = (lo - v) / d
            t2 = (hi - v) / d
        t_near = np.minimum(t1, t2)
        t_far = np.maximum(t1, t2)
        parallel = d == 0
        if np.any(parallel):
            inside = (v >= lo) & (v <= hi)
            t_near = np.where(parallel, np.where(inside, -np.inf, np.inf), t_near)
            t_far  = np.where(parallel, np.where(inside,  np.inf, -np.inf), t_far)
        return t_near, t_far

    tx_n, tx_f = slab_t(vx, dx, *box["x_range"])
    ty_n, ty_f = slab_t(vy, dy, *box["y_range"])
    tz_n, tz_f = slab_t(vz, dz, *box["z_range"])
    t_enter = np.maximum.reduce([tx_n, ty_n, tz_n, np.zeros_like(vx)])
    t_exit  = np.minimum.reduce([tx_f, ty_f, tz_f])
    return np.maximum(0.0, t_exit - t_enter)

DEFAULT_GSIMPLE_DIR = (
    "/cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr/sbnd/persistent/stash/"
    "fluxFiles/bnb/BooNEtoGSimple/configK-v1/july2023/neutrinoMode/"
)
DEFAULT_OUTPUT = "gsimple_flux.pkl"

PDG_BY_FLAVOR = {"nue": 12, "nuebar": -12, "numu": 14, "numubar": -14}
FLAVORS = list(PDG_BY_FLAVOR.keys())

E_BINS = np.linspace(0.0, 10.0, 81)
BIN_CENTERS = 0.5 * (E_BINS[:-1] + E_BINS[1:])


def _get_pot(f):
    meta_key = next((k for k in f.keys() if "meta" in k.lower()), None)
    if meta_key is None:
        raise RuntimeError("No meta tree found; cannot read POT.")
    meta = f[meta_key]
    pot_key = next((k for k in meta.keys() if "proton" in k.lower()), None)
    if pot_key is None:
        raise RuntimeError(f"No 'protons' branch in meta tree; keys={list(meta.keys())}")
    return float(meta[pot_key].array(library="np").sum())


def load_gsimple(filename):
    with uproot.open(filename) as f:
        pot = _get_pot(f)
        tree = f["flux"]
        prefix = "entry/" if "entry/pdg" in tree.keys() else ""
        branches = [prefix + b for b in ["pdg", "wgt", "vtxx", "vtxy", "vtxz", "px", "py", "pz", "E"]]
        arr = tree.arrays(branches, library="np")

    return {
        "pot": pot,
        "pdg": arr[prefix + "pdg"],
        "wgt": arr[prefix + "wgt"],
        "vtxx": arr[prefix + "vtxx"] * 100.0,
        "vtxy": arr[prefix + "vtxy"] * 100.0,
        "vtxz": arr[prefix + "vtxz"] * 100.0,
        "px": arr[prefix + "px"],
        "py": arr[prefix + "py"],
        "pz": arr[prefix + "pz"],
        "E":  arr[prefix + "E"],
    }


def path_length_volume(data, boxes):
    L = np.zeros_like(data["vtxx"])
    for box in boxes:
        L += path_length_box(
            data["vtxx"], data["vtxy"], data["vtxz"],
            data["dx"], data["dy"], data["dz"], box,
        )
    return L


def volume_cm3(boxes):
    v = 0.0
    for b in boxes:
        v += ((b["x_range"][1] - b["x_range"][0])
              * (b["y_range"][1] - b["y_range"][0])
              * (b["z_range"][1] - b["z_range"][0]))
    return v



def process_file(fpath):
    """Load one file and return POT-weighted histogram sums per volume and flavor."""
    try:
        chunk = load_gsimple(fpath)
    except Exception as e:
        print(f"  [skip] {os.path.basename(fpath)}: {e}", flush=True)
        return None

    p_mag = np.sqrt(chunk["px"]**2 + chunk["py"]**2 + chunk["pz"]**2)
    chunk["dx"] = chunk["px"] / p_mag
    chunk["dy"] = chunk["py"] / p_mag
    chunk["dz"] = chunk["pz"] / p_mag

    result = {"pot": chunk["pot"], "hists": {}}
    for vname, boxes in RAYTRACE_VOLUME_DEFS:
        L = path_length_volume(chunk, boxes)
        result["hists"][vname] = {}
        for flav in FLAVORS:
            pdg = PDG_BY_FLAVOR[flav]
            m = (chunk["pdg"] == pdg) & (L > 0)
            h, _ = np.histogram(chunk["E"][m], bins=E_BINS, weights=chunk["wgt"][m] * L[m])
            result["hists"][vname][flav] = h
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gsimple-dir", default=DEFAULT_GSIMPLE_DIR, help="Directory containing gsimple ROOT files.")
    ap.add_argument("--n-files", type=int, default=None, help="Max number of files to process (default: all).")
    ap.add_argument("--nproc", type=int, default=None, help="Number of parallel workers (default: 15).")
    ap.add_argument("--output", default=DEFAULT_OUTPUT, help="Output pickle path.")
    args = ap.parse_args()

    all_files = sorted(glob.glob(os.path.join(args.gsimple_dir, "*.root")))
    if not all_files:
        raise SystemExit(f"No ROOT files found in {args.gsimple_dir}")

    files = all_files[:args.n_files] if args.n_files is not None else all_files
    nproc = args.nproc or 15
    print(f"Processing {len(files)} file(s) with {nproc} workers...")

    # Accumulate histogram sums across files — never hold all raw arrays in memory
    total_pot = 0.0
    hist_sums = {vname: {flav: np.zeros(len(BIN_CENTERS)) for flav in FLAVORS}
                 for vname, _ in RAYTRACE_VOLUME_DEFS}

    with multiprocessing.Pool(processes=nproc) as pool:
        for result in tqdm(pool.imap_unordered(process_file, files), total=len(files), unit="file"):
            if result is None:
                continue
            total_pot += result["pot"]
            for vname in hist_sums:
                for flav in FLAVORS:
                    hist_sums[vname][flav] += result["hists"][vname][flav]

    if total_pot == 0.0:
        raise SystemExit("No files loaded successfully.")

    print(f"Total POT: {total_pot:.4g}")

    spectra = {}
    for vname, boxes in RAYTRACE_VOLUME_DEFS:
        V = volume_cm3(boxes)
        spectra[vname] = {
            flav: hist_sums[vname][flav] * 1.0e4 / (V * total_pot)
            for flav in FLAVORS
        }

    out = {
        "spectra":     spectra,
        "bin_edges":   E_BINS,
        "bin_centers": BIN_CENTERS,
        "total_pot":   total_pot,
        "flavors":     FLAVORS,
        "volumes":     dict(RAYTRACE_VOLUME_DEFS),
    }

    with open(args.output, "wb") as f:
        pickle.dump(out, f)
    print(f"Saved flux to {args.output}")

    for vname in spectra:
        nue_total = spectra[vname]["nue"] + spectra[vname]["nuebar"]
        print(f"  {vname}: nue+nuebar integrated flux = {nue_total.sum():.4e} /m^2/POT")


if __name__ == "__main__":
    main()
