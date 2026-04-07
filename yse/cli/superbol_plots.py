#!/usr/bin/env python3
"""Bolometric luminosity + interpolated multi-band plots (replaces five per-object superbol notebooks)."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from astropy.table import Table

from yse.paths import SUPERBOL_ROOT

# (subdir under superbol-master, output folder name, logL file, interp file, title)
OBJECTS = [
    ("20jfo", "superbol_output_2020jfo", "logL-obs_2020jfo_SDAuBgVriz.txt", "interpolated-lcs_2020jfo_SDAuBgVriz.txt", "2020jfo, ref = r"),
    ("20jww", "superbol_output_2020jww", "logL-obs_2020jww_gwrizy.txt", "interpolated-lcs_2020jww_gwrizy.txt", "2020jww Bol. Lum."),
    ("20rth", "superbol_output_2020rth", "logL-obs_2020rth_gcwroiz.txt", "interpolated-lcs_2020rth_gcwroiz.txt", "2020rth Bol. Lum."),
    ("20tly", "superbol_output_2020tly", "logL-obs_2020tly_gcwroiz.txt", "interpolated-lcs_2020tly_gcwroiz.txt", "2020tly Bol. Lum."),
    ("20hgw", "superbol_output_2020hgw", "logL-obs_2020hgw_SDAUBgVriz.txt", "interpolated-lcs_2020hgw_SDAUBgVriz.txt", "2020hgw Bol. Lum."),
]


def plot_one(tag: str, save_dir: Path | None) -> None:
    row = next((r for r in OBJECTS if r[0] == tag), None)
    if row is None:
        raise SystemExit(f"Unknown tag {tag}; choose one of: {[r[0] for r in OBJECTS]}")
    sub, outd, logf, intpf, ttl = row
    wdir = SUPERBOL_ROOT / sub / outd
    log_path = wdir / logf
    interp_path = wdir / intpf
    if not log_path.is_file():
        raise SystemExit(f"Missing {log_path}")

    table = Table.read(str(log_path), format="ascii", names=("MJD", "lum", "dlum"))
    plt.figure()
    plt.errorbar(table["MJD"], table["lum"], table["dlum"])
    plt.title(ttl)
    plt.xlabel("MJD")
    plt.ylabel("Luminosity")
    if save_dir:
        save_dir.mkdir(parents=True, exist_ok=True)
        p = save_dir / f"{tag}_logl.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        print(f"Wrote {p}")
    else:
        plt.show()
    plt.close()

    if not interp_path.is_file():
        print(f"Skip interp plot, missing {interp_path}")
        return

    interp_lc = Table.read(str(interp_path), format="ascii")
    colnames = list(interp_lc.colnames)
    phase_col = next(c for c in colnames if "phase" in str(c).lower())
    times = interp_lc[phase_col]
    rest = [c for c in colnames if c != phase_col]
    plt.figure(figsize=(10, 11))
    for i in range(0, len(rest), 2):
        band = rest[i]
        errcol = rest[i + 1] if i + 1 < len(rest) else None
        if band not in interp_lc.colnames:
            continue
        yerr = interp_lc[errcol] if errcol and errcol in interp_lc.colnames else None
        plt.errorbar(times, interp_lc[band], yerr=yerr, marker="o", label=str(band))
    plt.gca().invert_yaxis()
    plt.legend(ncol=5)
    plt.xlabel("MJD")
    plt.ylabel("App Mag")
    plt.title(ttl)
    if save_dir:
        p = save_dir / f"{tag}_interp_bands.png"
        plt.savefig(p, dpi=150, bbox_inches="tight")
        print(f"Wrote {p}")
    else:
        plt.show()
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("tag", nargs="?", help="Object folder tag, e.g. 20jfo; default: run all")
    p.add_argument("--save-dir", type=Path, help="Directory for PNG figures")
    args = p.parse_args()
    tags = [args.tag] if args.tag else [r[0] for r in OBJECTS]
    for t in tags:
        plot_one(t, args.save_dir)


if __name__ == "__main__":
    main()
