"""Shared formatting for extrabol ``.dat`` files and superbol per-filter text exports."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

# SVO-style filter strings for extrabol-style inputs (union of maps used across notebooks)
EXTRABOL_FLT_TO_SVO: dict[str, str] = {
    "B": "Swift/UVOT.B_fil",
    "U": "Swift/UVOT.U_fil",
    "UVM2": "Swift/UVOT.UVM2_fil",
    "UVW1": "Swift/UVOT.UVW1_fil",
    "UVW2": "Swift/UVOT.UVW2_fil",
    "U_S": "Swift/UVOT.U_fil",
    "V": "Swift/UVOT.V_fil",
    "G": "GAIA/GAIA0.G",
    "F555W": "HST/WFC3_UVIS1.F555W",
    "F814W": "HST/WFC3_UVIS1.F814W",
    "cyan-ATLAS": "Misc/Atlas.cyan",
    "orange-ATLAS": "Misc/Atlas.orange",
    "g": "PAN-STARRS/PS1.g",
    "r": "PAN-STARRS/PS1.r",
    "i": "PAN-STARRS/PS1.i",
    "w": "PAN-STARRS/PS1.w",
    "y": "PAN-STARRS/PS1.y",
    "z": "PAN-STARRS/PS1.z",
    "g-ZTF": "Palomar/ZTF.g",
    "r-ZTF": "Palomar/ZTF.r",
    "gp": "SLOAN/SDSS.gprime_filter",
    "rp": "SLOAN/SDSS.rprime_filter",
    "ip": "SLOAN/SDSS.iprime_filter",
    "up": "SLOAN/SDSS.uprime_filter",
    "i-N": "Generic/Cousins.I",
    "X": "Palomar/ZTF.g",
    "Y": "Palomar/ZTF.r",
}

# Single-letter / instrument tags expected by superbol file naming
SUPERBOL_FLT_REPLACEMENTS: dict[str, str] = {
    "g-ZTF": "g",
    "r-ZTF": "r",
    "orange-ATLAS": "o",
    "cyan-ATLAS": "c",
    "G": "E",
    "UVW1": "A",
    "UVM2": "D",
    "UVW2": "S",
    "i-N": "i",
    "gp": "g",
    "ip": "i",
    "rp": "r",
    "up": "u",
    "U_S": "u",
}


def extrabol_format_table(
    df: pd.DataFrame,
    *,
    flt_column: str = "FLT",
    flt_svo_map: dict[str, str] | None = None,
    magerr_max: float = 0.3,
    magerr_min: float = 0.0,
    drop_columns: Iterable[str] = ("FLUXCAL", "FLUXCALERR", "TELESCOPE", "MAGSYS"),
) -> pd.DataFrame:
    """Return MJD / MAG / MAGERR / FLT_SVO_ID / MAG TYPE columns for extrabol."""
    flt_svo_map = flt_svo_map or EXTRABOL_FLT_TO_SVO
    out = df.rename(columns={flt_column: "FLT_SVO_ID"}).copy()
    out["FLT_SVO_ID"] = out["FLT_SVO_ID"].astype(str).map(flt_svo_map)
    for c in drop_columns:
        if c in out.columns:
            out.drop(columns=[c], inplace=True)
    out["MAG TYPE"] = "AB"
    keep = [c for c in ("MJD", "MAG", "MAGERR", "FLT_SVO_ID", "MAG TYPE") if c in out.columns]
    out = out[keep].dropna(subset=["FLT_SVO_ID"])
    mask = (out["MAGERR"] < magerr_max) & (out["MAGERR"] > magerr_min)
    return out.loc[mask]


def write_extrabol_dat(df: pd.DataFrame, path: str | Path) -> Path:
    path = Path(path)
    df.to_csv(path, sep=" ", header=False, index=False)
    return path


def snana_passband_extrabol_table(
    df: pd.DataFrame,
    *,
    passband_to_svo: dict[str, str] | None = None,
    drop_phot: tuple[str, ...] = ("FLUX", "FLUXERR", "PHOTFLAG"),
) -> pd.DataFrame:
    """Light curves from :func:`yse.snana.read_YSE_ZTF_snana_dir` (PASSBAND column)."""
    pb = passband_to_svo or {
        "g": "PAN-STARRS/PS1.g",
        "r": "PAN-STARRS/PS1.r",
        "i": "PAN-STARRS/PS1.i",
        "z": "PAN-STARRS/PS1.z",
        "X": "Palomar/ZTF.g",
        "Y": "Palomar/ZTF.r",
    }
    out = df.rename(columns={"PASSBAND": "FLT_SVO_ID"}).copy()
    out["FLT_SVO_ID"] = out["FLT_SVO_ID"].map(pb)
    out = out.drop(columns=[c for c in drop_phot if c in out.columns])
    out["MAG TYPE"] = "AB"
    return out[["MJD", "MAG", "MAGERR", "FLT_SVO_ID", "MAG TYPE"]].dropna(subset=["FLT_SVO_ID"])


def write_superbol_filter_files(
    name: str,
    df: pd.DataFrame,
    *,
    out_dir: str | Path = ".",
    flt_column: str = "FLT",
    replacements: dict[str, str] | None = None,
    skip_filters: frozenset[str] | None = frozenset({"o", "c"}),
    magerr_max: float | None = 1.0,
) -> list[Path]:
    """
    Write ``{name}_{filter}.txt`` files with ``MJD MAG MAGERR`` per filter (superbol convention).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    replacements = replacements or SUPERBOL_FLT_REPLACEMENTS
    work = df.replace(replacements)
    cols = [c for c in (flt_column, "MJD", "MAG", "MAGERR") if c in work.columns]
    if flt_column not in cols:
        raise KeyError(f"Missing column {flt_column!r}")
    work = work[cols]
    if magerr_max is not None:
        work = work[work["MAGERR"] <= magerr_max]
    written: list[Path] = []
    for filt in work[flt_column].unique():
        fs = str(filt)
        if skip_filters and fs in skip_filters:
            continue
        part = work[work[flt_column] == filt][["MJD", "MAG", "MAGERR"]]
        path = out_dir / f"{name}_{fs}.txt"
        part.to_csv(path, sep=" ", header=False, index=False)
        written.append(path)
    return written
