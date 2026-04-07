#!/usr/bin/env python3
"""DR1 SNANA population utilities (replaces ``RSGlightcurves`` / ``typeiisorting`` notebooks)."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

# Editable install registers ``yse``; running from a clone adds the repo root.
_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from yse.paths import DR1_SNANA_II  # noqa: E402
from yse.snana import read_YSE_ZTF_snana_dir  # noqa: E402


def cmd_cuts(args: argparse.Namespace) -> None:
    snana_dir = Path(args.snana_dir).expanduser().resolve()
    full_snid_list, full_meta_list, full_df_list = read_YSE_ZTF_snana_dir(snana_dir, keep_ztf=True)
    first_cut = []
    for meta in full_meta_list:
        if meta["nobs_before_peak"] > 0 and meta["nobs_before_peak"] + meta["nobs_after_peak"] >= 7:
            first_cut.append(meta["object_id"])
    print(f"First cut (pre-peak + total length): {len(first_cut)} objects")

    all_df = []
    for df in full_df_list:
        u = df.replace(to_replace="X", value="g").replace(to_replace="Y", value="r")
        all_df.append(u)

    thesis_names: list[str] = []
    for df, name in zip(all_df, full_snid_list):
        peak = df["MAG"].min()
        peak_mjd = df.loc[df["MAG"] == peak, "MJD"].iloc[0]
        before_peak = df[df["MJD"] < peak_mjd]
        after_peak = df[df["MJD"] >= peak_mjd]
        ok = True
        for band in ("g", "r"):
            if before_peak[before_peak["PASSBAND"] == band].shape[0] < 5:
                ok = False
            if after_peak[after_peak["PASSBAND"] == band].shape[0] < 10:
                ok = False
        if ok:
            thesis_names.append(name)
    print(f"Thesis-style g/r depth cuts: {len(thesis_names)} objects")
    if args.print_names:
        for n in sorted(thesis_names):
            print(n)


def cmd_summary(args: argparse.Namespace) -> None:
    root = Path(args.snana_dir).expanduser().resolve()
    broad_ii: list[str] = []
    not_broad: list[str] = []
    subtype_counts: dict[str, int] = defaultdict(int)
    for path in sorted(root.glob("*.dat")):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        spec_broad = None
        for line in lines:
            if line.startswith("SPEC_CLASS_BROAD: "):
                spec_broad = line.split()
                break
        if spec_broad and len(spec_broad) > 2 and spec_broad[2] == "II":
            broad_ii.append(path.name)
        else:
            not_broad.append(path.name)
        if len(lines) > 19:
            subtype = lines[19].split()
            key = subtype[2] if len(subtype) > 2 else "?"
            subtype_counts[key] += 1
    print(f"SPEC_CLASS_BROAD == II: {len(broad_ii)} files")
    print(f"Other / imposters: {len(not_broad)} files")
    print("Subtype token (line ~20) histogram:")
    for k in sorted(subtype_counts):
        print(f"  {k}: {subtype_counts[k]}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("cuts", help="RSG / thesis-style light-curve cuts")
    pc.add_argument("snana_dir", type=Path, nargs="?", default=DR1_SNANA_II)
    pc.add_argument("--print-names", action="store_true")
    pc.set_defaults(func=cmd_cuts)

    ps = sub.add_parser("summary", help="Header / subtype histogram")
    ps.add_argument("snana_dir", type=Path, nargs="?", default=DR1_SNANA_II)
    ps.set_defaults(func=cmd_summary)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
