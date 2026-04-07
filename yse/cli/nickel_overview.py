#!/usr/bin/env python3
"""2020jfo multiband photometry plot + bolometric curve comparison (Nickel data tree)."""
# flake8: noqa

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from yse.paths import NICKEL_DIR, THESIS_PHOT
from yse.util import plot_multiband_lightcurve, read_whitespace_lightcurve


def plot_jfo_thesis_photometry() -> None:
    tab = read_whitespace_lightcurve(THESIS_PHOT / "2020jfo_lightcurves")
    df = tab.to_pandas()
    if "filter" in df.columns and "FLT" not in df.columns:
        df = df.rename(columns={"filter": "FLT"})
    plt.figure(figsize=(9, 6))
    plot_multiband_lightcurve(df, title="2020jfo (thesis photometry)", magerr_max=1.0)


def compare_bolometric_curves(
    ck_obs: Path,
    ck_bb: Path,
    extrabol_tab: Path,
    *,
    t0_mjd: float = 58975.25,
) -> None:
    ck = pd.read_csv(ck_obs, sep="\t", names=["MJD", "lum", "dlum"])
    ck_bb = pd.read_csv(ck_bb, sep="\t", names=["MJD", "lum", "dlum"])
    ex = pd.read_csv(extrabol_tab, sep=r"\s+")
    ex_t = ex["Time (MJD)"].astype(float)
    ex_l = np.log10(ex["Log10(Bol. Lum)"].astype(float))
    ex_e = ex["Log10(Bol. Err)"].astype(float) / ex["Log10(Bol. Lum)"].astype(float)
    plt.figure(figsize=(8, 6))
    plt.errorbar(ck["MJD"] - t0_mjd, ck["lum"], yerr=ck["dlum"], fmt="o", label="CK obs")
    plt.errorbar(ex_t, ex_l, yerr=ex_e, fmt="o", color="C3", label="extrabol")
    plt.xlabel("Phase (days)")
    plt.ylabel("log L")
    plt.legend()
    plt.figure(figsize=(8, 6))
    plt.errorbar(ck_bb["MJD"] - t0_mjd, ck_bb["lum"], yerr=ck_bb["dlum"], fmt="o", label="CK BB")
    plt.errorbar(ex_t, ex_l, yerr=ex_e, fmt="o", color="C3", label="extrabol")
    plt.xlabel("Phase (days)")
    plt.ylabel("log L")
    plt.legend()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phot", action="store_true", help="Plot thesis 2020jfo multiband")
    p.add_argument("--bol", action="store_true", help="Compare bolometric files under nickel data/yse")
    args = p.parse_args()
    yse_nickel = NICKEL_DIR / "data" / "yse"
    if args.phot:
        plot_jfo_thesis_photometry()
    if args.bol:
        compare_bolometric_curves(
            yse_nickel / "logL_obs_2020jfo_SDAUBgVriz.txt",
            yse_nickel / "logL_bb_2020jfo_SDAUBgVriz.txt",
            yse_nickel / "jfo_bol_LC",
        )
    if not args.phot and not args.bol:
        plot_jfo_thesis_photometry()
    plt.show()


if __name__ == "__main__":
    main()
