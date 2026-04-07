#!/usr/bin/env python3
"""
CLI for legacy bundled tools.

* ``smooth`` — Gaussian smoothing (replaces ``smoother_j2014.py`` / ``smoother_j2018.py``).

The Gaussian-process *extrabol* pipeline is :mod:`yse.extrabol.extrabol_gp` (run with
``python -m yse.extrabol.extrabol_gp``). Batch superbol plots are in
:mod:`yse.cli.superbol_plots`; upstream superbol inputs live under
``yse.paths.SUPERBOL_ROOT``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from yse.util.spectrum import smooth_flux_300kms


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    sm = sub.add_parser("smooth", help="Apply 300 km/s (default kernel) smoothing to a 2-column spectrum")
    sm.add_argument("path", type=Path, help="Input ASCII spectrum (wavelength, flux)")
    sm.add_argument("-o", "--output", type=Path, default=None, help="Output path (default: *_smoothed300kms.txt)")

    args = p.parse_args()
    if args.cmd == "smooth":
        out = smooth_flux_300kms(args.path, args.output)
        print(out)


if __name__ == "__main__":
    main()
