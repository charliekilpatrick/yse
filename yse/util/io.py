"""ASCII / Astropy table helpers for YSE-style photometry exports."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.table import Table


def read_yse_ascii_photometry(path: str | Path) -> Table:
    """
    Read YSE-PZ style ``*.snana.txt`` with a ``VARLIST:`` line and ``OBS:`` data rows.
    """
    path = Path(path)
    names: list[str] | None = None
    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("VARLIST:"):
            names = line.split()[1:]
        elif names is not None and line.startswith("OBS:"):
            parts = line.split()[1:]
            if len(parts) >= len(names):
                rows.append(parts[: len(names)])
            elif len(parts) == len(names) - 1:
                rows.append(parts)
    if not names or not rows:
        raise ValueError(f"No VARLIST/OBS table found in {path}")
    cols = list(zip(*rows))
    tab = Table([np.array(c) for c in cols], names=names)
    for col in ("MJD", "FLUXCAL", "FLUXCALERR", "MAG", "MAGERR"):
        if col in tab.colnames:
            tab[col] = np.array(tab[col], dtype=float)
    if "MAGERR" in tab.colnames:
        m = np.isfinite(tab["MAGERR"])
        tab = tab[m]
    return tab


def read_whitespace_lightcurve(path: str | Path) -> Table:
    """Load a space-separated photometry table with a header row (e.g. ``MJD FLT MAG MAGERR ...``)."""
    import pandas as pd

    df = pd.read_csv(path, sep=r"\s+", engine="python")
    return Table.from_pandas(df)


def list_snana_txt(directory: str | Path) -> list[Path]:
    d = Path(directory)
    return sorted(d.glob("*.snana.txt"))


def _snana_stem_key(path: Path) -> str:
    n = path.name
    if n.endswith("_data.snana.txt"):
        return n.replace("_data.snana.txt", "")
    if "_YSEdata.snana.txt" in n:
        return n.split("_")[0]
    return path.stem


def load_thesis_yse_photometry_tables(phot_dir: str | Path) -> dict[str, Table]:
    """
    All ``*.snana.txt`` tables under *phot_dir* plus ``2020jfo_lightcurves`` if present.
    Keys are object ids such as ``2020hgw``, ``2020tly``.
    """
    phot_dir = Path(phot_dir)
    out: dict[str, Table] = {}
    for p in list_snana_txt(phot_dir):
        out[_snana_stem_key(p)] = read_yse_ascii_photometry(p)
    jfo_path = phot_dir / "2020jfo_lightcurves"
    if jfo_path.is_file():
        out["2020jfo"] = read_whitespace_lightcurve(jfo_path)
    return out


def unique_filters(table: Table, column: str = "FLT") -> np.ndarray:
    return np.unique(np.array(table[column].astype(str)))
