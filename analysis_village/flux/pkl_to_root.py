#!/usr/bin/env python3
"""
Convert gsimple_flux.pkl to a ROOT file with per-flavor TH1D histograms
suitable for use with GENIE, NEUT, and NuWro (via PrepareGENIE / nuisflat).

Histogram names written:
  flux_sbnd_nue    (nue,     PDG  12)
  flux_sbnd_anue   (nuebar,  PDG -12)
  flux_sbnd_numu   (numu,    PDG  14)
  flux_sbnd_anumu  (numubar, PDG -14)

Units in the pickle: m^-2 POT^-1 per 50 MeV bin.
Units written to ROOT: cm^-2 POT^-1 per bin  (divide by 1e4).
The x-axis is neutrino energy in GeV.
"""

import argparse
import pickle

import numpy as np
import uproot

FLAVOR_TO_HIST = {
    "nue":     "flux_sbnd_nue",
    "nuebar":  "flux_sbnd_anue",
    "numu":    "flux_sbnd_numu",
    "numubar": "flux_sbnd_anumu",
}

M2_TO_CM2 = 1.0e-4


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input",  default="gsimple_flux.pkl",
                    help="Input pickle produced by compute_gsimple_flux.py")
    ap.add_argument("--volume", default="FV_split_truncY_eastonly",
                    help="Volume key to extract from the pickle")
    ap.add_argument("--output", default="sbnd_flux.root",
                    help="Output ROOT file path")
    args = ap.parse_args()

    with open(args.input, "rb") as f:
        data = pickle.load(f)

    if args.volume not in data["spectra"]:
        available = list(data["spectra"].keys())
        raise SystemExit(f"Volume '{args.volume}' not in pickle. Available: {available}")

    edges = data["bin_edges"].astype(np.float64)   # GeV, shape (81,)
    spectra = data["spectra"][args.volume]

    with uproot.recreate(args.output) as out:
        for flavor, hist_name in FLAVOR_TO_HIST.items():
            values = spectra[flavor].astype(np.float64) * M2_TO_CM2
            out[hist_name] = (values, edges)
            total = values.sum()
            print(f"  {hist_name:25s}  integral = {total:.4e} cm^-2 POT^-1")

    print(f"\nWrote {len(FLAVOR_TO_HIST)} histograms to {args.output}")
    print(f"Volume: {args.volume}  |  POT in sample: {data['total_pot']:.4g}")


if __name__ == "__main__":
    main()
