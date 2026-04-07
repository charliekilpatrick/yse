#!/usr/bin/env python3
"""Regenerate figures from the legacy ``nickel_analysis`` notebook.

Reads tables and BLAST JSON under ``yse.paths.NICKEL_DATA``, writes PNGs and
optional paired CSV summaries under ``--out-dir`` (default
``NICKEL_ANALYSIS_FIGURES``).

BLAST stage: every combination of SN quantities (Ni mass, peak luminosity, plateau
duration) × global/local host properties (27 panels, SN on *x*, host on *y*). After
all panels run, logs a rank-ordered list of Pearson *r* (by ``|r|``) and writes
``sn_host_correlation_rankings.txt`` beside the PNGs.

Example::

    python -m yse.cli.host_plots
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Mapping

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.cosmology import Planck15
from astropy.table import Column, Table, vstack
from scipy.stats import spearmanr

from yse.util.logging_config import configure_logging
from yse.paths import NICKEL_ANALYSIS_FIGURES, NICKEL_DATA

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Column / axis metadata (raw strings: LaTeX)
# ---------------------------------------------------------------------------

KEYWORDS: dict[str, dict[str, str]] = {
    "plateau_duration": {
        "key": "Tpt",
        "errkey": "e_Tpt",
        "latex": r"Rest-frame plateau duration (days)",
    },
    "nickel_mass": {
        "key": "MNi",
        "errkey": "e_MNi",
        "latex": r"Nickel-56 Ejecta Mass (M$_{\odot}$)",
    },
    "host_galaxy": {"key": "Host"},
    "peak_luminosity": {
        "key": "logLpeak",
        "errkey": "e_logLpeak",
        "latex": r"Peak Luminosity (log(erg s$^{-1}$))",
    },
}

BLAST_KEYWORDS: dict[str, dict[str, str]] = {
    "local_sfd": {
        "latex": r"$\Sigma$SFR ($M_\odot$ yr$^{-1}$ kpc$^{-2}$; log scale)",
    },
    "global_aperture_host_log_mass": {
        "latex": r"Global Stellar Mass (log($M_\odot$))",
    },
    "global_aperture_host_log_sfr": {
        "latex": r"Global Star-Formation Rate (log($M_\odot$ yr$^{-1}$))",
    },
    "global_aperture_host_log_ssfr": {
        "latex": r"Global Specific Star-Formation Rate (log(yr$^{-1}$))",
    },
    "local_aperture_host_logzsol": {
        "latex": r"Local Metallicity (log($Z/Z_\odot$))",
    },
    "global_aperture_host_logzsol": {
        "latex": r"Global Metallicity (log($Z/Z_\odot$))",
    },
    "local_aperture_host_log_age": {
        "latex": r"Local Host Stellar Age (Gyr)",
    },
    "local_aperture_host_logsfh": {
        "latex": r"Local Host Stellar Age (Gyr)",
    },
    "global_aperture_host_logsfh": {
        "latex": r"Global Host Stellar Age (Gyr)",
    },
}

COSMO = Planck15

STARRED_NAMES = ("2020hgw", "2020jfo", "2020jww", "1992af", "2002hx")

COLORS: dict[str, str] = {
    "valenti": "blue",
    "anderson": "blue",
    "yse": "blue",
    "peculiar": "red",
    "martinez": "blue",
    "spiro": "blue",
}

# Pearson correlation: parametric bootstrap (Gaussian measurement errors) for 95% CI
_CORR_BOOTSTRAP_N = 2000
_CORR_BOOTSTRAPRNG = np.random.default_rng(42)


def _pearson_r_xy(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson r on paired arrays (finite values only)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 2:
        return float("nan")
    if np.std(x, ddof=0) == 0.0 or np.std(y, ddof=0) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _spearman_r_xy(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation (robust to outliers); ignores NaNs in pairs."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r, _ = spearmanr(x, y, nan_policy="omit")
    return float(r) if np.isfinite(r) else float("nan")


def _pearson_correlation_with_errors(
    x: np.ndarray,
    y: np.ndarray,
    sx: np.ndarray,
    sy: np.ndarray,
    *,
    n_bootstrap: int = _CORR_BOOTSTRAP_N,
    rng: np.random.Generator | None = None,
) -> tuple[float, float, float]:
    """Pearson r on observed data and 95% CI from parametric error bootstrap.

    Each bootstrap draw adds independent Gaussian noise N(0, σ_x_i) and N(0, σ_y_i)
    to the measured points and recomputes Pearson r, so the interval reflects both
    sampling uncertainty and quoted measurement uncertainties (symmetric σ).

    Returns
    -------
    r_obs, ci_low, ci_high
    """
    rng = rng or _CORR_BOOTSTRAPRNG
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sx = np.asarray(sx, dtype=float)
    sy = np.asarray(sy, dtype=float)
    r_obs = _pearson_r_xy(x, y)
    n = x.size
    rs = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        xp = x + rng.normal(0.0, sx, size=n)
        yp = y + rng.normal(0.0, sy, size=n)
        rs[b] = _pearson_r_xy(xp, yp)
    ok = np.isfinite(rs)
    if not np.any(ok):
        return r_obs, float("nan"), float("nan")
    rs = rs[ok]
    ci_lo, ci_hi = np.percentile(rs, [2.5, 97.5])
    return r_obs, float(ci_lo), float(ci_hi)


def _correlation_panel_text(
    r_obs: float,
    ci_lo: float,
    ci_hi: float,
    rho: float,
) -> str:
    """One-line mathtext: Pearson *r* (optional bootstrap CI) | Spearman ρ."""
    sep = r"  |  "
    if np.isfinite(ci_lo) and np.isfinite(ci_hi):
        pearson = rf"$r = {r_obs:.2f}$ [{ci_lo:.2f}, {ci_hi:.2f}]"
    elif np.isfinite(r_obs):
        pearson = rf"$r = {r_obs:.2f}$"
    else:
        pearson = ""
    spearman = rf"$\rho = {rho:.2f}$" if np.isfinite(rho) else ""
    if pearson and spearman:
        return pearson + sep + spearman
    return pearson or spearman


# Shared typography for literature + BLAST figures
_AXIS_LABEL_FONTSIZE = 20
_TICK_LABEL_FONTSIZE = 18
_CORR_TITLE_FONTSIZE = 13
_LEGEND_FONTSIZE = 17


def _dedupe_legend(ax: plt.Axes) -> None:
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), fontsize=_LEGEND_FONTSIZE)


def _apply_nickel_axis_style(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", which="major", labelsize=_TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="minor", labelsize=max(_TICK_LABEL_FONTSIZE - 2, 12))


def _correlation_title_literature(r_obs: float, rho: float) -> str:
    """One-line title: Pearson *r* | Spearman ρ (no error bootstrap on *r*)."""
    sep = r"  |  "
    pearson = rf"$r = {r_obs:.2f}$" if np.isfinite(r_obs) else ""
    spearman = rf"$\rho = {rho:.2f}$" if np.isfinite(rho) else ""
    if pearson and spearman:
        return pearson + sep + spearman
    return pearson or spearman


# ---------------------------------------------------------------------------
# Data loading (notebook cells 3–11, 12)
# ---------------------------------------------------------------------------


def load_valenti(data_root: Path) -> Table:
    valenti_phys = Table.read(
        str(data_root / "valenti16" / "valenti_phys_bolometric.tsv"),
        format="ascii",
        delimiter="|",
    )
    valenti_lcpm = Table.read(
        str(data_root / "valenti16" / "valenti_phys_lcpm.tsv"),
        format="ascii",
        delimiter="|",
    )
    valenti_meta = Table.read(
        str(data_root / "valenti16" / "valenti_metadata.tsv"),
        format="ascii",
        delimiter="|",
    )
    for col in valenti_meta.keys():
        if col not in valenti_phys.keys() and col not in ("LC", "LC_param"):
            col_data = []
            for row in valenti_phys:
                mask = row["Name"] == valenti_meta["Name"]
                if len(valenti_meta[mask]) == 1:
                    col_data.append(valenti_meta[mask][0][col])
                else:
                    col_data.append(None)
            valenti_phys.add_column(Column(col_data, name=col))
    for col in valenti_lcpm.keys():
        if col not in valenti_phys.keys() and col not in ("LC", "LC_param"):
            col_data = []
            for row in valenti_phys:
                mask = row["Name"] == valenti_lcpm["Name"]
                if len(valenti_lcpm[mask]) == 1:
                    col_data.append(valenti_lcpm[mask][0][col])
                else:
                    col_data.append(None)
            valenti_phys.add_column(Column(col_data, name=col))
    for key in ("MNi", "e_MNi", "Tpt", "e_Tpt"):
        data = []
        for v in valenti_phys[key]:
            if v is None or (hasattr(v, "mask") and np.ma.is_masked(v)):
                data.append(np.nan)
            else:
                try:
                    data.append(float(v))
                except (TypeError, ValueError):
                    data.append(np.nan)
        valenti_phys[key] = data
    valenti_phys.meta["sample_name"] = "Normal"
    return valenti_phys


def load_spiro(data_root: Path) -> Table:
    spiro_data = Table.read(str(data_root / "spiro14" / "spiro_nickel"), format="ascii")
    spiro_data.meta["sample_name"] = "Normal"
    return spiro_data


def load_anderson(data_root: Path) -> Table:
    anderson_data = Table.read(
        str(data_root / "anderson14" / "anderson_fit_data.txt"), format="ascii"
    )
    mass: list[float] = []
    err_mass: list[tuple[float, float]] = []
    tpt: list[float] = []
    err_tpt: list[float] = []
    for row in anderson_data:
        dat = row["MNi"]
        if "$" in dat:
            val = float(dat.split("$")[0])
            err = dat.split("$")[1]
            err = err.replace("{", "").replace("}", "").replace("^", "")
            upper, lower = err.split("_")
            upper = abs(float(upper))
            lower = abs(float(lower))
            mass.append(val)
            err_mass.append((upper, lower))
        else:
            mass.append(np.nan)
            err_mass.append((np.nan, np.nan))
        dat = row["OPTd"]
        if "(" in dat:
            val, err = dat.split("(")
            tpt.append(float(val))
            err_tpt.append(float(err.replace(")", "")))
        else:
            tpt.append(np.nan)
            err_tpt.append(np.nan)
    all_lum: list[float] = []
    elum: list[float] = []
    for row in anderson_data:
        dat = row["M_max"]
        dat = dat.replace(")", "").replace("*", "")
        if "(" not in dat:
            all_lum.append(np.nan)
            elum.append(np.nan)
            continue
        val, err = dat.split("(")
        err = float(err) * 0.4
        lum = np.log10(3.839e33 * 10 ** (-0.4 * (float(val) - 4.74)))
        all_lum.append(lum)
        elum.append(err)
    anderson_data.add_column(Column(all_lum, name="logLpeak"))
    anderson_data.add_column(Column(elum, name="e_logLpeak"))
    anderson_data.rename_column("SN", "Name")
    anderson_data["MNi"] = mass
    anderson_data.add_column(Column(err_mass, name="e_MNi"))
    anderson_data.add_column(Column(tpt, name="Tpt"))
    anderson_data.add_column(Column(err_tpt, name="e_Tpt"))
    names: list[str] = []
    hosts: list[str] = []
    with open(data_root / "anderson14" / "anderson_sample.txt", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            names.append(parts[0])
            hosts.append(parts[1])
    names_a = np.array(names)
    hosts_a = np.array(hosts)
    a_hosts: list[str] = []
    for row in anderson_data:
        mask = names_a == row["Name"]
        if len(names_a[mask]) == 1:
            a_hosts.append(hosts_a[mask][0])
        else:
            a_hosts.append("--")
    anderson_data.add_column(Column(a_hosts, name="Host"))
    anderson_data.meta["sample_name"] = "Anderson+2014"
    return anderson_data


def load_martinez(data_root: Path) -> Table:
    martinez_data = Table.read(
        str(data_root / "martinez22" / "table_results.txt"),
        format="ascii",
        delimiter="&",
    )
    mass: list[float] = []
    err_mass: list[tuple[float, float]] = []
    for row in martinez_data:
        dat = row["MNi"]
        if "$" in dat:
            val = float(dat.split("$")[0])
            err = dat.split("$")[1]
            err = err.replace("{", "").replace("}", "").replace("^", "")
            upper, lower = err.split("_")
            mass.append(val)
            err_mass.append((abs(float(upper)), abs(float(lower))))
        else:
            mass.append(np.nan)
            err_mass.append((np.nan, np.nan))
    martinez_data["MNi"] = mass
    martinez_data.add_column(Column(err_mass, name="e_MNi"))
    martinez_data.meta["sample_name"] = "Martinez+2022"
    return martinez_data


def load_yse(data_root: Path) -> Table:
    yse_data = Table.read(str(data_root / "yse" / "yse_data.txt"), format="ascii")
    for key in ("e_Tpt", "e_MNi"):
        data = []
        for row in yse_data:
            data.append([float(r) for r in row[key].split(",")])
        yse_data[key] = data
    yse_data.meta["sample_name"] = "YSE"
    return yse_data


def split_starred_peculiars(all_data: dict[str, Table]) -> None:
    cols = [
        "Name",
        "MNi",
        "e_MNi",
        "Tpt",
        "e_Tpt",
        "logLpeak",
        "e_logLpeak",
        "FMT",
    ]
    star_list: list[Table] = []
    for key in list(all_data.keys()):
        table = all_data[key]
        table["FMT"] = "o"
        starred_mask = np.isin(table["Name"], STARRED_NAMES)
        table["FMT"][starred_mask] = "*"
        starred_table = table[starred_mask]
        if len(starred_table) > 0:
            starred_table["e_Tpt"] = [
                np.array([v, v]) if np.ndim(v) == 0 else v for v in starred_table["e_Tpt"]
            ]
            star_list.append(starred_table[cols])
        mask = ~np.isin(table["Name"], STARRED_NAMES)
        all_data[key] = table[mask]
    if star_list:
        starred = vstack(star_list, metadata_conflicts="silent")
        starred.meta["sample_name"] = "peculiar"
        all_data["peculiar"] = starred


def build_all_data(data_root: Path) -> dict[str, Table]:
    logger.info("Loading literature tables from %s", data_root)
    all_data: dict[str, Table] = {}
    all_data["valenti"] = load_valenti(data_root)
    logger.info("  valenti16: %d rows", len(all_data["valenti"]))
    all_data["spiro"] = load_spiro(data_root)
    logger.info("  spiro14: %d rows", len(all_data["spiro"]))
    all_data["anderson"] = load_anderson(data_root)
    logger.info("  anderson14: %d rows", len(all_data["anderson"]))
    all_data["martinez"] = load_martinez(data_root)
    logger.info("  martinez22: %d rows", len(all_data["martinez"]))
    all_data["yse"] = load_yse(data_root)
    logger.info("  yse: %d rows", len(all_data["yse"]))
    split_starred_peculiars(all_data)
    if "peculiar" in all_data:
        logger.info("  split peculiar subsample: %d rows", len(all_data["peculiar"]))
    logger.info(
        "Literature merge: %d tables (keys=%s)",
        len(all_data),
        ", ".join(sorted(all_data.keys())),
    )
    return all_data


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def make_plot(
    all_data: Mapping[str, Table],
    param1: str,
    param2: str,
) -> tuple[plt.Figure, plt.Axes]:
    logger.info("Building scatter: %s vs %s", param1, param2)
    fig, ax = plt.subplots(figsize=(8, 8))
    key1 = KEYWORDS[param1]["key"]
    key2 = KEYWORDS[param2]["key"]
    ekey1 = KEYWORDS[param1]["errkey"]
    ekey2 = KEYWORDS[param2]["errkey"]
    pooled_x: list[float] = []
    pooled_y: list[float] = []
    for key, table in all_data.items():
        if key1 not in table.keys() or key2 not in table.keys():
            continue
        if ekey1 not in table.keys() or ekey2 not in table.keys():
            continue
        x = np.asarray(table[key1].data, dtype=float)
        y = np.asarray(table[key2].data, dtype=float)
        xerr = np.asarray(table[ekey1].data)
        yerr = np.asarray(table[ekey2].data)
        mask = ~np.isnan(x) & ~np.isnan(y)
        x = x[mask]
        y = y[mask]
        n = mask.size
        if xerr.ndim == 2:
            if xerr.shape[0] == n:
                xerr = xerr[mask].T
            else:
                xerr = xerr[:, mask]
        else:
            xerr = xerr[mask]
        if yerr.ndim == 2:
            if yerr.shape[0] == n:
                yerr = yerr[mask].T
            else:
                yerr = yerr[:, mask]
        else:
            yerr = yerr[mask]
        if len(x) == 0:
            continue
        pooled_x.extend(x.tolist())
        pooled_y.extend(y.tolist())
        fmt = "*" if key == "peculiar" else "o"
        label = "Peculiar" if key == "peculiar" else "Normal"
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            label=label,
            fmt=fmt,
            markersize=10 if key == "peculiar" else 6,
            color=COLORS.get(key, "blue"),
        )
    ax.set_xlabel(KEYWORDS[param1]["latex"], fontsize=_AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(KEYWORDS[param2]["latex"], fontsize=_AXIS_LABEL_FONTSIZE)
    if len(pooled_x) >= 2:
        xa = np.asarray(pooled_x, dtype=float)
        ya = np.asarray(pooled_y, dtype=float)
        lit_title = _correlation_title_literature(_pearson_r_xy(xa, ya), _spearman_r_xy(xa, ya))
        if lit_title:
            ax.set_title(lit_title, fontsize=_CORR_TITLE_FONTSIZE)
    _apply_nickel_axis_style(ax)
    return fig, ax


def add_mesa_nickel_theory(ax: plt.Axes, theory_dir: Path) -> None:
    files = sorted(theory_dir.glob("mesa_nickel*.txt"))
    for file in files:
        filebase = file.name
        value = float(filebase.split("_")[-1].replace(".txt", "").replace("foe", ""))
        theory_table = Table.read(str(file), names=("MNi", "Tpt"), format="ascii")
        ax.plot(
            theory_table["Tpt"],
            theory_table["MNi"],
            linestyle="dashed",
            label=rf"MESA, $E={value}$ foe",
        )


def annotate_mesa_panel(ax: plt.Axes) -> None:
    ax.text(63, 0.045, "2020jfo")
    ax.text(89.5, 0.125, "2020hgw")
    ax.text(93, 0.11, "2020jww")
    ax.text(109, 0.07, "2020rth")
    ax.text(54, 0.079, "1992af")
    ax.text(67.3, 0.053, "2002hx")


def get_blast_data(
    name: str,
    blast_dir: Path,
    *,
    suppress_output: bool = False,
) -> dict[str, Any] | None:
    fname = blast_dir / f"{name}.json"
    if fname.is_file():
        with open(fname, encoding="utf-8") as f:
            return json.load(f)
    url = f"https://blast.scimma.org/api/transient/get/{name}"
    logger.debug("Fetching BLAST JSON from API: %s", name)
    r = requests.get(url, timeout=60)
    if r.status_code == 200:
        data = json.loads(r.content)
        blast_dir.mkdir(parents=True, exist_ok=True)
        with open(fname, "w", encoding="utf-8") as f:
            json.dump(data, f)
        logger.debug("Cached BLAST JSON to %s", fname)
        return data
    if not suppress_output:
        logger.warning("Could not get BLAST data for %s (HTTP %s)", name, r.status_code)
    else:
        logger.debug("Could not get BLAST data for %s (HTTP %s)", name, r.status_code)
    return None


def _local_aperture_area_kpc2(
    semi_major_arcsec: float,
    semi_minor_arcsec: float,
    angular_diameter_distance_kpc: float,
) -> float:
    """Ellipse area (kpc²) from angular semiaxes (arcsec) and angular-diameter distance (kpc)."""
    rad_per_arcsec = np.pi / (180.0 * 3600.0)
    r_a = float(semi_major_arcsec) * rad_per_arcsec * float(angular_diameter_distance_kpc)
    r_b = float(semi_minor_arcsec) * rad_per_arcsec * float(angular_diameter_distance_kpc)
    return float(np.pi * r_a * r_b)


def get_local_sfd(
    data: Mapping[str, Any],
    *,
    suppress_output: bool = False,
) -> tuple[float, tuple[float, float], str]:
    """Local ΣSFR (M☉ yr⁻¹ kpc⁻²) from BLAST log₁₀(SFR) posteriors and ellipse area.

    Central value is SFR₅₀ / area. Uncertainties are the distances from that value to
    the 16th and 84th ΣSFR percentiles (linear SFR ÷ same area), so the band
    [Σ₁₆, Σ₈₄] is always non-negative. If posteriors are numerically non-monotonic,
    the three Σ values are sorted before taking differences.
    """
    error = ""
    redshift = data["host_redshift"]
    if redshift is None:
        name = data["transient_name"]
        if not suppress_output:
            logger.warning("Redshift missing for %s (local SFD)", name)
        else:
            logger.debug("Redshift missing for %s (local SFD)", name)
        return np.nan, (np.nan, np.nan), "redshift"
    a = data["local_aperture_semi_major_axis_arcsec"]
    b = data["local_aperture_semi_minor_axis_arcsec"]
    if a is None or b is None:
        name = data["transient_name"]
        if not suppress_output:
            logger.warning("No local aperture for %s (local SFD)", name)
        else:
            logger.debug("No local aperture for %s (local SFD)", name)
        return np.nan, (np.nan, np.nan), "local_aperture"

    d_kpc = float(COSMO.angular_diameter_distance(redshift).value) * 1000.0
    area_kpc2 = _local_aperture_area_kpc2(float(a), float(b), d_kpc)
    if not np.isfinite(area_kpc2) or area_kpc2 <= 0.0:
        return np.nan, (np.nan, np.nan), "area"

    log16 = data.get("local_aperture_host_log_sfr_16")
    log50 = data.get("local_aperture_host_log_sfr_50")
    log84 = data.get("local_aperture_host_log_sfr_84")
    if log16 is None or log50 is None or log84 is None:
        return np.nan, (np.nan, np.nan), "local_aperture"

    l16, l50, l84 = float(log16), float(log50), float(log84)
    if not all(np.isfinite(v) for v in (l16, l50, l84)):
        return np.nan, (np.nan, np.nan), "local_aperture"

    # Linear SFR (M☉ yr⁻¹) from log₁₀ posteriors, then ΣSFR = SFR / area (kpc²)
    sfr_16 = 10.0**l16
    sfr_50 = 10.0**l50
    sfr_84 = 10.0**l84
    sigma_16 = sfr_16 / area_kpc2
    sigma_50 = sfr_50 / area_kpc2
    sigma_84 = sfr_84 / area_kpc2

    if not (sigma_16 <= sigma_50 <= sigma_84):
        sigma_16, sigma_50, sigma_84 = sorted((sigma_16, sigma_50, sigma_84))
        name = data.get("transient_name", "?")
        if not suppress_output:
            logger.warning(
                "BLAST local log SFR percentiles not monotonic for %s; "
                "using sorted ΣSFR bounds for error bars",
                name,
            )
        else:
            logger.debug("Non-monotonic ΣSFR percentiles for %s; sorted bounds", name)

    err_upper = max(0.0, sigma_84 - sigma_50)
    err_lower = max(0.0, sigma_50 - sigma_16)
    return sigma_50, (err_upper, err_lower), error


def make_blast_plot(
    all_data: Mapping[str, Table],
    param1: str,
    param2: str,
    blast_dir: Path,
    *,
    suppress_output: bool = False,
    csv_path: Path | None = None,
) -> tuple[plt.Figure, plt.Axes, dict[str, str], float]:
    """Return figure, axes, BLAST error map, and Pearson *r* on pooled finite pairs (or NaN)."""
    fig, ax = plt.subplots(figsize=(8, 8))
    r_pearson_out = float("nan")
    blast1 = param1.startswith("global") or param1.startswith("local")
    blast2 = param2.startswith("global") or param2.startswith("local")
    blast_errors: dict[str, str] = {}

    if param1 in KEYWORDS:
        key1 = KEYWORDS[param1]["key"]
        ekey1 = KEYWORDS[param1]["errkey"]
        latex1 = KEYWORDS[param1]["latex"]
    else:
        latex1 = BLAST_KEYWORDS[param1]["latex"]
    if param2 in KEYWORDS:
        key2 = KEYWORDS[param2]["key"]
        ekey2 = KEYWORDS[param2]["errkey"]
        latex2 = KEYWORDS[param2]["latex"]
    else:
        latex2 = BLAST_KEYWORDS[param2]["latex"]

    all_x: list[float] = []
    all_y: list[float] = []
    all_xerr_upper: list[float] = []
    all_xerr_lower: list[float] = []
    all_yerr_upper: list[float] = []
    all_yerr_lower: list[float] = []

    for key, table in all_data.items():
        if not blast1 and key1 in table.keys() and ekey1 in table.keys():
            x = np.array(table[key1].data)
            xerr = np.array(table[ekey1].data)
        elif blast1:
            x_list: list[float] = []
            xerr_list: list[tuple[float, float]] = []
            for row in table:
                data = get_blast_data(row["Name"], blast_dir, suppress_output=suppress_output)
                if data is None:
                    blast_errors.setdefault(row["Name"], "missing")
                    x_list.append(np.nan)
                    xerr_list.append((np.nan, np.nan))
                    continue
                if param1 == "local_sfd":
                    val, unc, error = get_local_sfd(data, suppress_output=suppress_output)
                    if error:
                        blast_errors.setdefault(row["Name"], error)
                else:
                    pk = param1 + "_50"
                    if data.get(pk) is None:
                        x_list.append(np.nan)
                        xerr_list.append((np.nan, np.nan))
                        continue
                    val = data[pk]
                    eup = data[param1 + "_84"] - val
                    elo = val - data[param1 + "_16"]
                    unc = (eup, elo)
                x_list.append(float(val))
                # (lower, upper) for matplotlib yerr/xerr shape (2, N): row 0 = lower, row 1 = upper
                xerr_list.append((float(unc[1]), float(unc[0])))
            x = np.array(x_list, dtype=float)
            xerr = np.array(xerr_list, dtype=float)
        else:
            continue

        if not blast2 and key2 in table.keys() and ekey2 in table.keys():
            y = np.array(table[key2].data)
            yerr = np.array(table[ekey2].data)
        elif blast2:
            y_list: list[float] = []
            yerr_list: list[tuple[float, float]] = []
            for row in table:
                data = get_blast_data(row["Name"], blast_dir, suppress_output=suppress_output)
                if data is None:
                    blast_errors.setdefault(row["Name"], "missing")
                    y_list.append(np.nan)
                    yerr_list.append((np.nan, np.nan))
                    continue
                if param2 == "local_sfd":
                    val, unc, error = get_local_sfd(data, suppress_output=suppress_output)
                    if error:
                        blast_errors.setdefault(row["Name"], error)
                elif param2 == "local_aperture_host_logsfh":
                    block = data.get("local_aperture_host_logsfh")
                    if block is None:
                        y_list.append(np.nan)
                        yerr_list.append((np.nan, np.nan))
                        continue
                    b0 = block[0]
                    val = b0["logsfr_50"]
                    eup = b0["logsfr_84"] - val
                    elo = val - b0["logsfr_16"]
                    unc = (eup, elo)
                elif param2 == "global_aperture_host_logsfh":
                    block = data.get("global_aperture_host_logsfh")
                    if block is None:
                        y_list.append(np.nan)
                        yerr_list.append((np.nan, np.nan))
                        continue
                    b0 = block[0]
                    val = b0["logsfr_50"]
                    eup = b0["logsfr_84"] - val
                    elo = val - b0["logsfr_16"]
                    unc = (eup, elo)
                else:
                    pk = param2 + "_50"
                    if data.get(pk) is None:
                        y_list.append(np.nan)
                        yerr_list.append((np.nan, np.nan))
                        continue
                    val = data[pk]
                    eup = data[param2 + "_84"] - val
                    elo = val - data[param2 + "_16"]
                    unc = (eup, elo)
                y_list.append(float(val))
                yerr_list.append((float(unc[1]), float(unc[0])))
            y = np.array(y_list, dtype=float)
            yerr = np.array(yerr_list, dtype=float)
        else:
            continue

        xerr = np.array(xerr)
        yerr = np.array(yerr)
        if xerr.ndim == 2:
            xerr = xerr.transpose()
        if yerr.ndim == 2:
            yerr = yerr.transpose()

        x1 = np.array(x, dtype=float)
        y1 = np.array(y, dtype=float)
        x1err = np.array(xerr)
        y1err = np.array(yerr)
        mask = ~np.isnan(x1) & ~np.isnan(y1)
        x1 = x1[mask]
        y1 = y1[mask]

        all_x.extend(x1.tolist())
        all_y.extend(y1.tolist())

        if x1err.ndim == 2:
            x1err = x1err[:, mask]
            all_xerr_lower.extend(x1err[0, :].tolist())
            all_xerr_upper.extend(x1err[1, :].tolist())
        else:
            x1err = x1err[mask]
            all_xerr_upper.extend(x1err.tolist())
            all_xerr_lower.extend(x1err.tolist())

        if y1err.ndim == 2:
            y1err = y1err[:, mask]
            all_yerr_lower.extend(y1err[0, :].tolist())
            all_yerr_upper.extend(y1err[1, :].tolist())
        else:
            y1err = y1err[mask]
            all_yerr_upper.extend(y1err.tolist())
            all_yerr_lower.extend(y1err.tolist())

        fmt = "*" if key == "peculiar" else "o"
        label = "Peculiar" if key == "peculiar" else "Normal"
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            label=label,
            fmt=fmt,
            markersize=10 if key == "peculiar" else 6,
            color=COLORS.get(key, "blue"),
        )

    all_x_arr = np.array(all_x, dtype=float)
    all_y_arr = np.array(all_y, dtype=float)
    all_xerr_upper = np.array(all_xerr_upper, dtype=float)
    all_xerr_lower = np.array(all_xerr_lower, dtype=float)
    all_yerr_upper = np.array(all_yerr_upper, dtype=float)
    all_yerr_lower = np.array(all_yerr_lower, dtype=float)

    if csv_path is not None and len(all_x_arr):
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write(
                f"{param1},{param1}_err_upper,{param1}_err_lower,"
                f"{param2},{param2}_err_upper,{param2}_err_lower\n"
            )
            for i in range(len(all_x_arr)):
                f.write(
                    f"{all_x_arr[i]},{all_xerr_upper[i]},{all_xerr_lower[i]},"
                    f"{all_y_arr[i]},{all_yerr_upper[i]},{all_yerr_lower[i]}\n"
                )
        logger.info("Wrote ODR sample table %s (%d points)", csv_path, len(all_x_arr))

    if blast_errors:
        logger.info(
            "BLAST panel %s vs %s: %d objects with missing/partial BLAST data",
            param1,
            param2,
            len(blast_errors),
        )

    if len(all_x_arr) >= 2:
        all_xerr_val = np.array(
            [max(all_xerr_upper[i], all_xerr_lower[i]) for i in range(len(all_x_arr))]
        )
        all_yerr_val = np.array(
            [max(all_yerr_upper[i], all_yerr_lower[i]) for i in range(len(all_y_arr))]
        )

        # Pearson r + 95% CI (parametric bootstrap, Gaussian σ_x/σ_y); Spearman ρ on ranks.
        eps = np.finfo(float).eps
        sx = np.maximum(all_xerr_val, eps * 1e6)
        sy = np.maximum(all_yerr_val, eps * 1e6)
        r_obs, ci_lo, ci_hi = _pearson_correlation_with_errors(
            all_x_arr, all_y_arr, sx, sy
        )
        r_pearson_out = r_obs
        rho = _spearman_r_xy(all_x_arr, all_y_arr)
        panel_txt = _correlation_panel_text(r_obs, ci_lo, ci_hi, rho)
        if panel_txt:
            ax.set_title(panel_txt, fontsize=_CORR_TITLE_FONTSIZE)
        logger.debug(
            "Correlation %s vs %s: r=%s [%s, %s], rho=%s (n=%d)",
            param1,
            param2,
            f"{r_obs:.4f}" if np.isfinite(r_obs) else "nan",
            f"{ci_lo:.4f}" if np.isfinite(ci_lo) else "nan",
            f"{ci_hi:.4f}" if np.isfinite(ci_hi) else "nan",
            f"{rho:.4f}" if np.isfinite(rho) else "nan",
            len(all_x_arr),
        )

    ax.set_xlabel(latex1, fontsize=_AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(latex2, fontsize=_AXIS_LABEL_FONTSIZE)
    if param1 == "local_sfd":
        ax.set_xscale("log", nonpositive="clip")
    if param2 == "local_sfd":
        ax.set_yscale("log", nonpositive="clip")
    _apply_nickel_axis_style(ax)
    return fig, ax, blast_errors, r_pearson_out


# ---------------------------------------------------------------------------
# SN (literature) × BLAST host: every combination (x = SN quantity, y = host quantity)
# ---------------------------------------------------------------------------

SN_HOST_PARAMETERS: tuple[str, ...] = (
    "nickel_mass",
    "peak_luminosity",
    "plateau_duration",
)

HOST_BLAST_PARAMETERS: tuple[str, ...] = (
    "global_aperture_host_log_mass",
    "global_aperture_host_log_sfr",
    "global_aperture_host_log_ssfr",
    "global_aperture_host_logzsol",
    "global_aperture_host_logsfh",
    "local_sfd",
    "local_aperture_host_logzsol",
    "local_aperture_host_log_age",
    "local_aperture_host_logsfh",
)

BLAST_PANELS: list[tuple[str, str]] = [
    (sn, host) for sn in SN_HOST_PARAMETERS for host in HOST_BLAST_PARAMETERS
]


def _sn_host_rank_file_entry(
    rank: int,
    stem: str,
    sn_key: str,
    host_key: str,
    r_pearson: float,
) -> str:
    sn_lab = KEYWORDS[sn_key]["latex"].replace("$", "")
    host_lab = BLAST_KEYWORDS[host_key]["latex"].replace("$", "")
    return (
        f"{rank:3d}.  r={r_pearson:+.4f}  |r|={abs(r_pearson):.4f}  file={stem}.png\n"
        f"      SN parameter ({sn_key}, x-axis): {sn_lab}\n"
        f"      Host parameter ({host_key}, y-axis): {host_lab}"
    )


def write_sn_host_correlation_rankings(
    out_dir: Path,
    entries: list[tuple[str, str, str, float]],
) -> None:
    """Log and save Pearson *r* rank (by |r|) for SN vs host panels.

    *entries* — ``(stem, sn_key, host_key, r_pearson)`` per panel.
    """
    finite = [(s, sn, h, r) for s, sn, h, r in entries if np.isfinite(r)]
    finite.sort(key=lambda t: (-abs(t[3]), t[0]))
    lines: list[str] = [
        "SN–host Pearson correlation rank (strongest |r| first). "
        "x-axis = SN quantity; y-axis = BLAST host quantity.",
        "",
    ]
    logger.info(
        "SN–host correlation ranking by |Pearson r| (%d of %d panels with finite r):",
        len(finite),
        len(entries),
    )
    for i, (stem, sn_key, host_key, r_pearson) in enumerate(finite, start=1):
        logger.info(
            "%3d. r=%+.4f  |r|=%.4f  %s vs %s",
            i,
            r_pearson,
            abs(r_pearson),
            sn_key,
            host_key,
        )
        lines.append(_sn_host_rank_file_entry(i, stem, sn_key, host_key, r_pearson))
        lines.append("")
    out_path = out_dir / "sn_host_correlation_rankings.txt"
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info("Wrote ranking list to %s", out_path)


def savefig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("Wrote %s", path)


def run_all_figures(
    out_dir: Path,
    data_root: Path | None = None,
    *,
    blast_only: bool = False,
    literature_only: bool = False,
    write_csv: bool = True,
    suppress_blast_log: bool = True,
) -> None:
    root = data_root or NICKEL_DATA
    blast_dir = root / "blast"
    theory_dir = root / "theory"
    logger.info(
        "Starting host_plots (data_root=%s, out_dir=%s, blast_only=%s, "
        "literature_only=%s, write_csv=%s)",
        root,
        out_dir,
        blast_only,
        literature_only,
        write_csv,
    )
    all_data = build_all_data(root)
    png_dir = out_dir
    csv_dir = out_dir / "tables"

    if not blast_only:
        logger.info("Literature figures: peak/MNi and plateau/MNi+MESA")
        fig, ax = make_plot(all_data, "peak_luminosity", "nickel_mass")
        _dedupe_legend(ax)
        savefig(fig, png_dir / "peak_luminosity_vs_nickel_mass.png")

        fig, ax = make_plot(all_data, "plateau_duration", "nickel_mass")
        n_theory = len(list(theory_dir.glob("mesa_nickel*.txt")))
        logger.info("Adding MESA theory tracks from %s (%d files)", theory_dir, n_theory)
        add_mesa_nickel_theory(ax, theory_dir)
        annotate_mesa_panel(ax)
        _dedupe_legend(ax)
        savefig(fig, png_dir / "plateau_duration_vs_nickel_mass_mesa.png")

    if not literature_only:
        n_panels = len(BLAST_PANELS)
        logger.info("BLAST panels: %d SN×host plots (cache/API under %s)", n_panels, blast_dir)
        rank_entries: list[tuple[str, str, str, float]] = []
        for i, (param1, param2) in enumerate(BLAST_PANELS, start=1):
            stem = f"{param1}_vs_{param2}"
            logger.info("BLAST panel %d/%d: %s", i, n_panels, stem)
            csv_path = csv_dir / f"{stem}.csv" if write_csv else None
            fig, ax, _, r_pearson = make_blast_plot(
                all_data,
                param1,
                param2,
                blast_dir,
                suppress_output=suppress_blast_log,
                csv_path=csv_path,
            )
            rank_entries.append((stem, param1, param2, r_pearson))
            _dedupe_legend(ax)
            savefig(fig, png_dir / f"{stem}.png")
        write_sn_host_correlation_rankings(png_dir, rank_entries)

    logger.info("Finished host_plots (output under %s)", out_dir)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out-dir",
        type=Path,
        default=NICKEL_ANALYSIS_FIGURES,
        help=f"PNG (and optional CSV) output directory (default: {NICKEL_ANALYSIS_FIGURES})",
    )
    p.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Override nickel data root (default: NICKEL_DATA)",
    )
    p.add_argument(
        "--literature",
        action="store_true",
        help="Only peak/MNi and plateau/MNi+MESA (default: full suite)",
    )
    p.add_argument(
        "--blast",
        action="store_true",
        help="Only BLAST correlation panels (default: full suite)",
    )
    p.add_argument("--no-csv", action="store_true", help="Skip writing ODR sample tables")
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Log DEBUG messages (API/cache detail) to stdout",
    )
    p.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Only log WARNING and above (less progress on stdout)",
    )
    p.add_argument(
        "--verbose-blast",
        action="store_true",
        help="Promote BLAST/local-SFR issues to WARNING (stderr); default is DEBUG only",
    )
    args = p.parse_args(argv)
    if args.verbose and args.quiet:
        p.error("Choose at most one of --verbose and --quiet")
    log_level = logging.WARNING if args.quiet else logging.DEBUG if args.verbose else logging.INFO
    configure_logging(level=log_level)

    literature = args.literature
    blast = args.blast
    if literature and blast:
        p.error("Choose at most one of --literature and --blast")
    blast_only = blast and not literature
    literature_only = literature and not blast
    if not literature and not blast:
        literature_only = blast_only = False

    run_all_figures(
        args.out_dir,
        args.data_root,
        blast_only=blast_only,
        literature_only=literature_only,
        write_csv=not args.no_csv,
        suppress_blast_log=not args.verbose_blast,
    )


if __name__ == "__main__":
    main()
