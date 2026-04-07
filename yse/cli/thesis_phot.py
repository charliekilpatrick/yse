#!/usr/bin/env python3
"""Thesis ``all_YSE_phot``: inspect tables, plot, write extrabol/superbol inputs."""
# flake8: noqa

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt

from yse import THESIS_PHOT
from yse.util import (
    extrabol_format_table,
    list_snana_txt,
    load_thesis_yse_photometry_tables,
    plot_multiband_lightcurve,
    read_yse_ascii_photometry,
    unique_filters,
    write_extrabol_dat,
    write_superbol_filter_files,
)

DEFAULT_OBJECTS = ("2020hgw", "2020jfo", "2020jww", "2020rth", "2020tly")


def cmd_check(_: argparse.Namespace) -> None:
    tly = read_yse_ascii_photometry(THESIS_PHOT / "2020tly_YSEdata.snana.txt")
    print("columns:", tly.colnames)
    print("filters:", unique_filters(tly))


def cmd_list(_: argparse.Namespace) -> None:
    for f in list_snana_txt(THESIS_PHOT):
        print(f.name)


def cmd_plot(args: argparse.Namespace) -> None:
    tables = load_thesis_yse_photometry_tables(THESIS_PHOT)
    names = args.objects or DEFAULT_OBJECTS
    for i, name in enumerate(names):
        df = tables[name].to_pandas().dropna()
        plt.figure(i)
        plot_multiband_lightcurve(df, title=name, magerr_max=args.magerr_max)
    plt.show()


def cmd_extrabol(args: argparse.Namespace) -> None:
    tables = load_thesis_yse_photometry_tables(THESIS_PHOT)
    names = args.objects or DEFAULT_OBJECTS
    for name in names:
        df = tables[name].to_pandas().dropna()
        out = extrabol_format_table(df, flt_column="FLT")
        write_extrabol_dat(out, f"{name}.dat")
        print(name, len(out), "rows")


def cmd_superbol(args: argparse.Namespace) -> None:
    tables = load_thesis_yse_photometry_tables(THESIS_PHOT)
    names = args.objects or DEFAULT_OBJECTS
    for name in names:
        df = tables[name].to_pandas().dropna()
        df = df[df["MAGERR"] <= args.magerr_max]
        n = len(write_superbol_filter_files(name, df, magerr_max=None, skip_filters=None))
        print(name, n, "filter files")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("check", help="Print columns/filters for 2020tly").set_defaults(func=cmd_check)
    sub.add_parser("list", help="List *.snana.txt in thesis photometry dir").set_defaults(func=cmd_list)

    pp = sub.add_parser("plot", help="Multiband figures for named objects")
    pp.add_argument(
        "--objects",
        nargs="*",
        default=list(DEFAULT_OBJECTS),
        help=f"Object ids (default: {DEFAULT_OBJECTS})",
    )
    pp.add_argument("--magerr-max", type=float, default=0.3, dest="magerr_max")
    pp.set_defaults(func=cmd_plot)

    pe = sub.add_parser("extrabol", help="Write extrabol-style .dat files")
    pe.add_argument("--objects", nargs="*", default=list(DEFAULT_OBJECTS))
    pe.set_defaults(func=cmd_extrabol)

    ps = sub.add_parser("superbol", help="Write per-filter superbol text files")
    ps.add_argument("--objects", nargs="*", default=list(DEFAULT_OBJECTS))
    ps.add_argument("--magerr-max", type=float, default=1.0, dest="magerr_max")
    ps.set_defaults(func=cmd_superbol)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
