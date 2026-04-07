#!/usr/bin/env python3
"""
Valenti-style bolometric fit + emcee + 56Ni scaling vs SN 1987A for thesis SNe.
Replaces the near-duplicate ``rth_nickel`` / ``jww_nickel`` / ``hgw_nickel`` notebooks.
"""

from __future__ import annotations

import argparse
import copy
from pathlib import Path

import emcee
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Column, Table

from yse.paths import NICKEL_DIR

# Bolometric light-curve files live under ``nickel/data/yse/`` in the dataset.
_BOL = NICKEL_DIR / "data" / "yse"

# (data file name under nickel/data/yse/, mask_date, ref_mjd_for_ni, plot title, print_extra_chi)
CASES = {
    "rth": ("rth_bol_LC", 15, 110.0, "20rth Valenti fit", False),
    "jww": ("jww_bol_LC", 20, 105.0, "20jww Valenti fit", False),
    "hgw": ("hgw_bol_LC", 20, 110.0, "20hgw Valenti fit", True),
}


def valenti_model(time: np.ndarray, theta) -> np.ndarray:
    a, tpt, w, p, m = theta
    return -a / (1.0 + np.exp((time - tpt) / w)) + p * time + m


def lnlikelihood(theta, lum, dlum, time) -> float:
    model_lum = valenti_model(time, theta)
    a, tpt, w, p, m = theta
    if a < 0.01 or a > 10.0 or tpt < 10.0 or w < 0.01:
        return -np.inf
    inv_sigma = 1.0 / dlum**2
    return -0.5 * np.sum((lum - model_lum) ** 2 * inv_sigma)


def load_table(path: Path, mjd_start: float, mask_date: float) -> Table:
    table = Table.read(str(path), format="ascii", names=("MJD", "lum", "dlum"))
    mbol = -2.5 * np.log10(table["lum"]) + 83.9605452803 + 4.74
    mbolerr = 1.086 * (table["dlum"] / table["lum"])
    table.add_column(Column(mbol, name="Mbol"))
    table.add_column(Column(mbolerr, name="Mbolerr"))
    table.add_column(Column(table["MJD"].data - mjd_start, name="t_exp"))
    table = table[table["t_exp"] > mask_date]
    print(table)
    return table


def mcmc_valenti(path: Path, mjd_start: float, mask_date: float, inflate_errors: float = 0.1):
    ndim = 5
    nwalkers = 100
    params = [1.7391707763335567, 111.06179709805089, 7.203197652533321, 0.008194307783125725, -13.583601384015793]
    pos = [params + 1e-2 * np.random.randn(ndim) for _ in range(nwalkers)]
    table = load_table(path, mjd_start, mask_date)
    fit_table = table.copy()
    fit_table = fit_table[np.abs(fit_table["MJD"] - mjd_start) >= 1]
    fit_table["Mbolerr"] = np.sqrt(fit_table["Mbolerr"] ** 2 + inflate_errors**2)
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        lnlikelihood,
        args=(fit_table["Mbol"], fit_table["Mbolerr"], fit_table["t_exp"]),
    )
    sampler.run_mcmc(pos, 5000, progress=True)
    return sampler.chain[:, 99:, :].reshape((-1, ndim)), table


def get_param_uncertainties(samples) -> list:
    vals = [
        v
        for v in map(
            lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
            zip(*np.percentile(samples, [16, 50, 84], axis=0)),
        )
    ]
    return vals


def erg(x: float) -> float:
    lsun = 3.838e33
    return lsun * 10 ** (0.4 * (4.74 - x))


def abs_mag_from_lum(tbl: Table) -> np.ndarray:
    mbol_sun = 4.74
    log_l_sun = 83.9605452803
    return -2.5 * tbl["lum"] + log_l_sun + mbol_sun


def run_case(case: str, mjd_start: float, savefig: Path | None) -> None:
    stem, mask_date, ref_day, title, extra_chi = CASES[case]
    data_path = _BOL / stem
    if not data_path.is_file():
        raise SystemExit(f"Missing {data_path}")

    samples, table = mcmc_valenti(data_path, mjd_start=mjd_start, mask_date=mask_date)
    params = get_param_uncertainties(samples)
    print("Best-fit parameters:")
    for p in copy.copy(params):
        print(f"Param: {p[0]} + {p[1]} - {p[2]}")
    best_fit = [p[0] for p in copy.copy(params)]
    print(best_fit)
    plt.figure()
    plt.errorbar(
        table["MJD"] - mjd_start,
        table["Mbol"],
        yerr=table["Mbolerr"],
        marker="*",
        zorder=5,
        label="data",
    )
    model_lum = valenti_model(table["MJD"] - mjd_start, best_fit)
    plt.plot(table["MJD"] - mjd_start, model_lum, zorder=10, label="valenti model")
    plt.ylabel("Absolute Magnitude")
    plt.xlabel("Time from Explosion")
    plt.legend()
    plt.title(title)
    plt.gca().invert_yaxis()
    if savefig:
        plt.savefig(savefig, dpi=150, bbox_inches="tight")
        print(f"Wrote {savefig}")
    else:
        plt.show()

    lightcurve_ref_time = ref_day - mjd_start
    print(lightcurve_ref_time)
    luminosity = []
    for s in samples:
        abs_mag = valenti_model(lightcurve_ref_time, s)
        luminosity.append(3.839e33 * 10 ** (-0.4 * (abs_mag - 4.74)))
    lum_pct = np.percentile(luminosity, [16, 50, 84])
    print("Luminosity percentiles:", lum_pct)

    epsilon = 6.8e9
    mass_co56 = np.array(luminosity) / epsilon
    lambda_co56 = np.log(2) / 55.9383
    lambda_ni56 = np.log(2) / 6.075
    mass_ni56 = mass_co56 * (lambda_co56 - lambda_ni56) / lambda_ni56 * (
        np.exp(-lambda_ni56 * lightcurve_ref_time) - np.exp(-lambda_co56 * lightcurve_ref_time)
    ) ** -1
    m_pct = np.percentile(mass_ni56, [16, 50, 84])
    print("56Ni mass percentiles (g):", m_pct)

    if extra_chi:
        chi_sq = np.sum((model_lum - table["Mbol"]) ** 2 / table["Mbolerr"] ** 2)
        print("chi_sq", chi_sq, "Nobs", len(table))

    chisq = float(np.sum((model_lum - table["Mbol"]) ** 2 / table["Mbolerr"] ** 2))
    dof = len(table["Mbol"]) - 5
    print("reduced chi^2", chisq / dof)

    sn87a_path = NICKEL_DIR / "data" / "1987A.dat"
    sn87a = Table.read(str(sn87a_path), format="ascii", names=("MJD", "lum", "dlum"))
    mbol_87a = abs_mag_from_lum(sn87a)
    sn_interp = np.interp(ref_day, table["MJD"], table["Mbol"])
    l_sn = erg(sn_interp)
    sn87a_interp = np.interp(ref_day, sn87a["MJD"], mbol_87a)
    l_87a = erg(float(sn87a_interp))
    ni_87a = 0.075
    ni_est = ni_87a * (l_sn / l_87a)
    print(f"Ni mass scaled from 87A (~Msun): {ni_est / 1.989e33:.4f}")

    peak_l = float(np.max(table["lum"]))
    peak_idx = int(np.argmax(table["lum"]))
    peak_err = (peak_l - table["dlum"][peak_idx]) / peak_l
    print(f"Peak L = {peak_l} +/- {peak_err}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case", choices=sorted(CASES), help="Which SN bolometric series to fit")
    p.add_argument("--mjd-start", type=float, default=0.0, help="Explosion/reference MJD offset")
    p.add_argument("--savefig", type=Path, help="Save fit figure instead of showing")
    args = p.parse_args()
    run_case(args.case, args.mjd_start, args.savefig)


if __name__ == "__main__":
    main()
