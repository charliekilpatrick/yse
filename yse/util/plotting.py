"""Matplotlib helpers for multiband photometry (thesis-style SNANA tables)."""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Common passband → color map for mixed-instrument thesis photometry
THESIS_BAND_COLORS: dict[str, str] = {
    "B": "blue",
    "F555W": "yellowgreen",
    "F814W": "darkred",
    "G": "black",
    "U": "purple",
    "UVM2": "darkviolet",
    "UVW1": "magenta",
    "UVW2": "plum",
    "U_S": "orchid",
    "V": "silver",
    "cyan-ATLAS": "cyan",
    "g-ZTF": "green",
    "g": "olive",
    "gp": "darkolivegreen",
    "i": "yellow",
    "i-N": "gold",
    "ip": "lightyellow",
    "orange-ATLAS": "orange",
    "r": "red",
    "r-ZTF": "crimson",
    "rp": "firebrick",
    "up": "fuchsia",
    "w": "teal",
    "y": "navy",
    "z": "lightsteelblue",
    "X": "blue",
    "Y": "yellow",
}


def plot_multiband_lightcurve(
    df: pd.DataFrame,
    *,
    flt_col: str = "FLT",
    mjd_col: str = "MJD",
    mag_col: str = "MAG",
    magerr_col: str = "MAGERR",
    colors: dict[str, str] | None = None,
    magerr_max: float | None = 0.3,
    require_nonneg_magerr: bool = True,
    title: str = "",
    ax: Any | None = None,
    invert_y: bool = True,
    legend_ncol: int = 4,
) -> Any:
    """Plot every distinct filter in *df* as errorbar points on one axes."""
    colors = colors or THESIS_BAND_COLORS
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 10))
    for f in sorted(df[flt_col].astype(str).unique()):
        m = df[flt_col].astype(str) == f
        if magerr_max is not None:
            m &= df[magerr_col] < magerr_max
        if require_nonneg_magerr:
            m &= df[magerr_col] >= 0
        sub = df.loc[m]
        if sub.empty:
            continue
        c = colors.get(str(f), "#666666")
        ax.errorbar(
            sub[mjd_col],
            sub[mag_col],
            yerr=sub[magerr_col],
            fmt="o",
            ms=4,
            color=c,
            markeredgecolor="k",
            label=str(f),
        )
    if invert_y:
        ax.invert_yaxis()
    ax.set_xlabel("Modified Julian Date")
    ax.set_ylabel("Brightness (mag)")
    if title:
        ax.set_title(title)
    ax.legend(ncol=legend_ncol, fontsize="small")
    return ax


def plot_dr1_snana_bands(
    df: pd.DataFrame,
    title: str,
    *,
    passband_col: str = "PASSBAND",
    xy_filters: tuple[str, ...] = ("X", "Y"),
    optical_filters: tuple[str, ...] = ("g", "r", "i", "z"),
    magerr_max: float = 0.2,
) -> Any:
    """DR1 ``*.dat`` dataframe with PAN-STARRS + survey X/Y columns."""
    _, ax = plt.subplots()
    cmap = {"X": "blue", "Y": "yellow", "g": "green", "r": "red", "i": "orange", "z": "purple"}
    for filt in xy_filters + optical_filters:
        m = (df[passband_col] == filt) & (df["MAGERR"] < magerr_max)
        sub = df.loc[m]
        if sub.empty:
            continue
        ax.errorbar(
            sub["MJD"],
            sub["MAG"],
            yerr=sub["MAGERR"],
            fmt="o",
            color=cmap.get(filt, "0.5"),
            markeredgecolor="k",
            label=filt,
        )
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel("MJD")
    ax.set_ylabel("AB mag")
    ax.legend()
    return ax


def plot_normalized_spectra_csv(
    paths: list,
    *,
    redshift: float,
    wavelength_col: str = "wavelength",
    flux_col: str = "flux",
    ascii_kwargs: dict | None = None,
) -> Any:
    """Overplot normalized spectra read with :func:`astropy.io.ascii.read`."""
    from astropy.io import ascii

    ascii_kwargs = ascii_kwargs or {"format": "csv", "header_start": 18, "data_start": 19}
    _, ax = plt.subplots()
    for p in paths:
        t = ascii.read(str(p), **ascii_kwargs)
        w = np.asarray(t[wavelength_col], dtype=float)
        fl = np.asarray(t[flux_col], dtype=float)
        med = np.median(fl) if np.any(fl) else 1.0
        ax.plot(w / (1.0 + redshift), fl / med, label=p.name if hasattr(p, "name") else str(p))
    ax.set_xlabel(r"Rest $\lambda$ (Å)")
    ax.set_ylabel("Normalized flux")
    ax.legend(fontsize=8)
    return ax
