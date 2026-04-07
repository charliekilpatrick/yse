#!/usr/bin/env python3
"""DR1 SNANA directory exploration and YSE-PZ-style griz depth selection on ascii photometry."""
# flake8: noqa

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt

from yse import DR1_SNANA_II, YSEPZ, read_YSE_ZTF_snana_dir
from yse.util import plot_dr1_snana_bands, select_snana_by_griz_depth, snana_passband_extrabol_table


def cmd_dr1(args: argparse.Namespace) -> None:
    root = Path(args.snana_dir or os.environ.get("YSE_SNANA_DIR", DR1_SNANA_II))
    snids, metas, dfs = read_YSE_ZTF_snana_dir(root)
    for m in metas:
        print(
            f"{m['object_id']} RA={m['ra']} DEC={m['dec']} z={m['redshift']} "
            f"host={m.get('host_gal_name')} gap={m.get('max_mjd_gap')}"
        )
    if args.plot:
        for snid, df in zip(snids, dfs):
            plt.figure()
            plot_dr1_snana_bands(df, snid)
        plt.show()
    if args.extrabol_preview:
        for snid, df in zip(snids, dfs):
            print(snid)
            print(snana_passband_extrabol_table(df).head())


def cmd_ysepz(args: argparse.Namespace) -> None:
    phot = Path(args.phot_dir or YSEPZ)
    kept, names = select_snana_by_griz_depth(
        phot,
        min_per_filter_before_peak=args.before,
        min_per_filter_after_peak=args.after,
    )
    print(f"kept {len(kept)} / files in {phot}")
    for n in names:
        print(" ", n)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("dr1", help="Load *.dat SNANA dir, print meta, optional plots")
    d.add_argument("--snana-dir", type=Path, default=None)
    d.add_argument("--plot", action="store_true")
    d.add_argument("--extrabol-preview", action="store_true")
    d.set_defaults(func=cmd_dr1)

    y = sub.add_parser("ysepz", help="Select objects in phot_dir passing griz pre/post peak cuts")
    y.add_argument("--phot-dir", type=Path, default=None, help=f"default: {YSEPZ}")
    y.add_argument("--before", type=int, default=5)
    y.add_argument("--after", type=int, default=10)
    y.set_defaults(func=cmd_ysepz)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
