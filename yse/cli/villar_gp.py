#!/usr/bin/env python3
"""II-P light curves: Villar/Bazin GP fit per band (needs ``light-curve`` + ``george``)."""
# flake8: noqa

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import astropy.table as at
try:
    import light_curve as lcpy
except ImportError as e:
    raise SystemExit("pip install light-curve george") from e

FILT_COL_MAP = {"g": "C2", "r": "C3", "i": "C1", "c": "cyan", "o": "gold", "z": "tab:grey", "w": "black"}
SURVEY_MARKERS = {"PAN-STARRS/PS1": "D", "Palomar/ZTF": "*", "Misc/Atlas": "v"}


def load_extrabol_style_lcs(pattern: str = "20*dat", cwd: Path | None = None) -> dict[str, at.Table]:
    cwd = cwd or Path.cwd()
    out: dict[str, at.Table] = {}
    for file in sorted(cwd.glob(pattern)):
        obj = file.name.replace(".dat", "").split("_")[0]
        lc = at.Table.read(
            file,
            names=("mjd", "mag", "mag_err", "pb", "zpt"),
            format="ascii.no_header",
            data_start=2,
        )
        lc = lc[lc["mag_err"] < 0.3]
        survey, filt = zip(*[x.split(".") for x in lc["pb"]])
        lc["survey"] = survey
        lc["filt"] = [x.replace("prime_filter", "").replace("cyan", "c").replace("orange", "o") for x in filt]
        out[obj] = lc
    return out


def fit_and_plot_object(obj: str, lc: at.Table, outdir: Path | None = None) -> None:
    outdir = outdir or Path.cwd()
    pb_counts = Counter(lc["filt"])
    fig, ax = plt.subplots(figsize=(8, 8))
    ctr = 0.0
    for pb in FILT_COL_MAP:
        this_pb = lc["filt"] == pb
        x = np.asarray(lc["mjd"][this_pb], dtype=float)
        m = np.asarray(lc["mag"][this_pb], dtype=float)
        merr = np.asarray(lc["mag_err"][this_pb], dtype=float)
        survey = lc["survey"][this_pb]
        if pb_counts[pb] <= 0:
            continue
        flag = 0
        for s in SURVEY_MARKERS:
            ind = survey == s
            if len(x[ind]) == 0:
                continue
            flag += 1
            label = f"{pb}+{ctr:0.1f}" if flag == 1 else None
            ax.errorbar(
                x[ind],
                m[ind] + ctr,
                yerr=merr[ind],
                marker=SURVEY_MARKERS[s],
                linestyle="none",
                label=label,
                color=FILT_COL_MAP[pb],
                ms=5,
            )
        if pb_counts[pb] <= 8:
            ctr += 0.5
            continue
        f = 10.0 ** (-0.4 * (m - 25.0))
        ferr = f * (np.log(10.0) / 2.5) * merr
        peakguess = float(x[np.argmax(f)])
        sel = (x > peakguess - 100.0) & (x < peakguess + 500.0)
        frange = float(f[sel].max() - f[sel].min())
        fmin = float(f[sel].min())
        rf = (f - fmin) / frange
        rferr = ferr / frange
        x, m, merr, rf, rferr, survey = x[sel], m[sel], merr[sel], rf[sel], rferr[sel], survey[sel]
        x_pred = np.arange(np.floor(x.min()), np.ceil(x.max()), 1.0)
        bazin = lcpy.VillarFit("mcmc", mcmc_niter=5000)
        bazinres = bazin(x, rf, sigma=rferr, sorted=True)
        bazinpred = bazin.model(x_pred, bazinres)
        predmag = -2.5 * np.log10((bazinpred * frange) + fmin) + 25.0
        ax.plot(x_pred, predmag + ctr, "-", color=FILT_COL_MAP[pb], alpha=0.5)
        ctr += 0.5
        ax.set_xlim(float(x.min()) - 5, float(x.max()) + 5)

    ax.legend(loc="upper right", frameon=False, fontsize="medium")
    ax.set_xlabel("MJD")
    ax.set_ylabel("mag + offset")
    ax.invert_yaxis()
    fig.suptitle(obj)
    fig.tight_layout()
    fig.savefig(outdir / f"{obj}_Villar_Fit.pdf")


def main() -> None:
    data = load_extrabol_style_lcs()
    for obj, lc in data.items():
        fit_and_plot_object(obj, lc)


if __name__ == "__main__":
    main()
