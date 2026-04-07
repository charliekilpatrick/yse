"""Sample cuts on light-curve tables (e.g. YSE-PZ style griz depth)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from yse.util.io import list_snana_txt, read_yse_ascii_photometry


def peak_mjd_from_min_mag(df: pd.DataFrame, mag_col: str = "MAG", mjd_col: str = "MJD") -> float:
    peak_mag = df[mag_col].min()
    return float(df.loc[df[mag_col] == peak_mag, mjd_col].iloc[0])


def griz_pre_post_peak_counts(
    df: pd.DataFrame,
    *,
    flt_col: str = "FLT",
    mjd_col: str = "MJD",
    mag_col: str = "MAG",
    filters: tuple[str, ...] = ("g", "r", "i", "z"),
) -> tuple[dict[str, int], dict[str, int]]:
    """Return counts per filter before / on-or-after peak MJD (griz unified names)."""
    work = df.replace({"g-ZTF": "g", "r-ZTF": "r"}).dropna()
    pk = peak_mjd_from_min_mag(work, mag_col=mag_col, mjd_col=mjd_col)
    before = work[work[mjd_col] < pk]
    after = work[work[mjd_col] >= pk]
    before_c = {f: int((before[flt_col] == f).sum()) for f in filters}
    after_c = {f: int((after[flt_col] == f).sum()) for f in filters}
    return before_c, after_c


def select_snana_by_griz_depth(
    phot_dir: str | Path,
    *,
    min_per_filter_before_peak: int = 5,
    min_per_filter_after_peak: int = 10,
    flt_col: str = "FLT",
) -> tuple[list[pd.DataFrame], list[str]]:
    """
    Load every ``*.snana.txt`` in *phot_dir* and keep objects meeting griz counts
    before/after peak (same logic as the old ysepz sample-selection notebook).
    """
    phot_dir = Path(phot_dir)
    names: list[str] = []
    kept: list[pd.DataFrame] = []
    for p in list_snana_txt(phot_dir):
        tab = read_yse_ascii_photometry(p)
        df = tab.to_pandas()
        b, a = griz_pre_post_peak_counts(df, flt_col=flt_col)
        if all(b[f] >= min_per_filter_before_peak for f in b) and all(
            a[f] >= min_per_filter_after_peak for f in a
        ):
            kept.append(df)
            names.append(p.name)
    return kept, names
