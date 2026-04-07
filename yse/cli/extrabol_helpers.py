#!/usr/bin/env python3
"""Small utilities: bolometric ascii tables, extrabol filter checks, CURVE peak-L helpers."""
# flake8: noqa

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.table import Column, Table


def extrabol_dat_unique_filters(path: Path, *, data_start: int = 2) -> np.ndarray:
    """Unique filter column (col4) in a 2-line-header extrabol ``.dat`` file."""
    t = Table.read(str(path), format="ascii", data_start=data_start)
    return np.unique(np.asarray(t["col4"]))


def load_bolometric_table(path: str | Path, mjd_start: float, mask_t_exp_gt: float) -> Table:
    """MJD / lum / dlum ascii → Mbol, t_exp columns (from old CURVE_peak_L notebook)."""
    path = Path(path)
    table = Table.read(path, format="ascii", names=("MJD", "lum", "dlum"))
    mbol = -2.5 * np.log10(table["lum"]) + 83.9605452803 + 4.74
    mbolerr = 1.086 * (table["dlum"] / table["lum"])
    table.add_column(Column(mbol, name="Mbol"))
    table.add_column(Column(mbolerr, name="Mbolerr"))
    t_exp = table["MJD"].data - mjd_start
    table.add_column(Column(t_exp, name="t_exp"))
    return table[table["t_exp"] > mask_t_exp_gt]


def plot_peak_luminosity_sampling(path: str | Path, *, stride: int = 5) -> None:
    t = Table.read(path, format="ascii", names=("MJD", "lum", "dlum"))
    lum = np.asarray(t["lum"], dtype=float)
    dlum = np.asarray(t["dlum"], dtype=float)
    yerr = (lum - dlum) / np.where(lum != 0, lum, 1.0)
    plt.figure()
    plt.errorbar(t["MJD"][::stride], lum[::stride], yerr=yerr[::stride])
    plt.ylabel("Lum")
    plt.xlabel("MJD")


def strip_extrabol_columns_to_bol_lc(in_path: str | Path) -> Path:
    """Drop temperature/radius columns from extrabol output; write ``*_bol_LC.txt``."""
    in_path = Path(in_path)
    bol = pd.read_csv(in_path, sep=r"\s+")
    drop = [c for c in bol.columns if "Temp" in c or "Radius" in c]
    bol = bol.drop(columns=drop, errors="ignore")
    stem = in_path.name.split("_")[0]
    out = in_path.with_name(stem + "_bol_LC.txt")
    bol.to_csv(out, sep=" ", index=False)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("filt-check", help="Print unique filters in extrabol .dat")
    f.add_argument("dat", type=Path)

    pk = sub.add_parser("peak-l-plot", help="Plot lum vs MJD from 3-column ascii")
    pk.add_argument("path", type=Path)
    pk.add_argument("--stride", type=int, default=5)

    args = p.parse_args()
    if args.cmd == "filt-check":
        print(extrabol_dat_unique_filters(args.dat))
    elif args.cmd == "peak-l-plot":
        plot_peak_luminosity_sampling(args.path, stride=args.stride)
        plt.show()


if __name__ == "__main__":
    main()
