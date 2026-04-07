"""Passband name normalization (e.g. YSE X/Y → SDSS-style g/r labels)."""

from __future__ import annotations

import pandas as pd


def uniformize_yse_xy_passbands(
    df: pd.DataFrame,
    column: str = "PASSBAND",
    x_to: str = "g",
    y_to: str = "r",
) -> pd.DataFrame:
    """Return a copy with ``X`` → *x_to* and ``Y`` → *y_to* in *column* (in-place safe)."""
    out = df.copy()
    if column not in out.columns:
        return out
    out[column] = out[column].replace({"X": x_to, "Y": y_to})
    return out
