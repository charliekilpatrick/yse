#!/usr/bin/env python3
"""Plot normalized spectra from CSVs in a working directory (Kast / LRIS style exports)."""
# flake8: noqa

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii

# Default redshifts for thesis sample objects (rest-frame plots)
REDSHIFTS = {"hgw": 0.043, "jfo": 0.005, "jww": 0.046, "rth": 0.0184, "tly": 0.058}


def plot_csv_glob(
    cwd: Path,
    *,
    pattern: str = "*.csv",
    z_default: float = 0.0,
    skip_substrings: tuple[str, ...] = ("jfo", "hgw"),
    ascii_kwargs: dict | None = None,
) -> None:
    ascii_kwargs = ascii_kwargs or {"format": "csv", "header_start": 18, "data_start": 19}
    for path in sorted(cwd.glob(pattern)):
        name = path.name.lower()
        if any(s in name for s in skip_substrings):
            continue
        z = z_default
        for key, zi in REDSHIFTS.items():
            if key in name:
                z = zi
                break
        t = ascii.read(str(path), **ascii_kwargs)
        w = np.asarray(t["wavelength"], dtype=float)
        fl = np.asarray(t["flux"], dtype=float)
        med = np.median(fl) if np.any(fl) else 1.0
        plt.figure()
        plt.plot(w / (1.0 + z), fl / med)
        plt.title(path.name)
        plt.xlabel(r"Rest $\lambda$ (Å)")
        plt.ylabel("Normalized flux")


def plot_na_d_region(
    path: Path,
    z: float,
    *,
    rest_line_aa: float = 5892.0,
    half_width_aa: float = 20.0,
    poly_deg: int = 2,
) -> None:
    """Local continuum + normalized flux around Na I D (one spectrum)."""
    t = ascii.read(str(path), format="csv", header_start=18, data_start=19)
    wave = np.asarray(t["wavelength"], dtype=float)
    flux = np.asarray(t["flux"], dtype=float)
    cen = rest_line_aa * (1.0 + z)
    lo, hi = cen - half_width_aa, cen + half_width_aa
    cont_lo = (wave > lo - 150) & (wave < lo - 10)
    cont_hi = (wave > hi + 10) & (wave < hi + 150)
    continuum = cont_lo | cont_hi
    coef = np.polyfit(wave[continuum], flux[continuum], poly_deg)
    cont = np.poly1d(coef)(wave)
    plt.figure()
    plt.plot(wave, flux / cont)
    plt.axvline(cen, color="0.5", ls="--")
    plt.xlim(lo - 80, hi + 80)
    plt.title(path.name)
    plt.xlabel(r"$\lambda$ (Å)")
    plt.ylabel("Normalized flux")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--cwd",
        type=Path,
        default=Path(os.environ.get("YSE_SPECTRA_CWD", ".")),
        help="Directory containing spectrum CSVs",
    )
    p.add_argument("--pattern", default="*.csv")
    p.add_argument("--na-demo", action="store_true", help="Plot Na D for one hard-coded path if present")
    args = p.parse_args()
    if args.na_demo:
        demo = args.cwd / "2020jww-KAST-2020-05-23.csv"
        if demo.is_file():
            plot_na_d_region(demo, REDSHIFTS["jww"])
        else:
            print("missing", demo)
    else:
        plot_csv_glob(args.cwd, pattern=args.pattern)
    plt.show()


if __name__ == "__main__":
    main()
