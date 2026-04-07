"""
Accelerated kernels for :mod:`yse.extrabol.extrabol_gp`.

Uses NumPy vectorization everywhere; optional Numba JIT (``pip install numba``)
speeds tight loops. Falls back cleanly when Numba is absent.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit, prange

    _HAS_NUMBA = True
except ImportError:
    _HAS_NUMBA = False

    def njit(*args, **kwargs):
        def deco(f):
            return f

        if args and callable(args[0]):
            return args[0]
        return deco

    def prange(x):
        return range(x)


# Physical constants (cgs), duplicated for Numba kernels
_H = 6.62607e-27
_C = 2.99792458e10
_K_B = 1.38064852e-16
_PI = np.pi
_ANG_TO_CM = 1e-8
_TWO_PI_H_C2 = 2.0 * _PI * _H * _C**2
_HC_KB = _H * _C / _K_B


@njit(cache=True, parallel=True)
def bbody_lam_flux_njit(lam_angstrom: np.ndarray, T: float, R: float) -> np.ndarray:
    """Blackbody L_lambda (erg/s/cm) at wavelengths ``lam_angstrom`` (Å); scalar T (K), R (cm)."""
    n = lam_angstrom.shape[0]
    out = np.empty(n, dtype=np.float64)
    area = 4.0 * _PI * R * R
    for i in prange(n):
        lam_cm = lam_angstrom[i] * _ANG_TO_CM
        expn = _HC_KB / (lam_cm * T)
        # guard overflow for very small T or large lam
        if expn > 700.0:
            blam = 0.0
        else:
            ex = np.exp(expn) - 1.0
            if ex <= 0.0:
                blam = 0.0
            else:
                blam = _TWO_PI_H_C2 / (lam_cm**5) / ex
        out[i] = blam * area
    return out


def bbody_numpy(lam, T, R):
    """Vectorized Planck L_lambda; same convention as legacy ``bbody``."""
    lam = np.asarray(lam, dtype=np.float64)
    scalar_in = lam.ndim == 0
    if scalar_in:
        lam = lam.reshape(1)
    if _HAS_NUMBA:
        out = bbody_lam_flux_njit(lam.ravel(), float(T), float(R))
        out = out.reshape(lam.shape)
    else:
        lam_cm = lam * _ANG_TO_CM
        exponential = (_H * _C) / (lam_cm * _K_B * T)
        blam = (_TWO_PI_H_C2 / (lam_cm**5)) / (np.exp(exponential) - 1.0)
        area = 4.0 * _PI * R**2
        out = blam * area
    if scalar_in:
        return float(out[0])
    return out


def chi_square_vec(dat: np.ndarray, model: np.ndarray, uncertainty: np.ndarray) -> float:
    """Sum of squared normalized residuals (replaces Python for-loop)."""
    u = np.asarray(uncertainty, dtype=np.float64)
    d = np.asarray(dat, dtype=np.float64)
    m = np.asarray(model, dtype=np.float64)
    return float(np.sum(((m - d) / u) ** 2))


def get_flam_batch(
    mag: np.ndarray,
    magerr: np.ndarray,
    wavelength_angstrom: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorized ``get_flam`` for many points (mag, magerr, wavelength in Å).
    Matches single-point ``get_flam`` in extrabol_gp.
    """
    mag = np.asarray(mag, dtype=np.float64)
    magerr = np.asarray(magerr, dtype=np.float64)
    wavelength = np.asarray(wavelength_angstrom, dtype=np.float64)
    fnu = 10.0 ** ((-mag + 48.6) / -2.5)
    fnu = fnu * 4.0 * np.pi * (3.086e19) ** 2
    fnu_err = (
        np.abs(0.921034 * 10.0 ** (0.4 * mag - 19.44))
        * magerr
        * 4.0
        * np.pi
        * (3.086e19) ** 2
    )
    wcm = wavelength * _ANG_TO_CM
    flam = fnu * _C / (wcm**2)
    flam_err = fnu_err * _C / (wcm**2)
    return flam, flam_err


def build_filter_index_map(filter_ids: np.ndarray) -> dict:
    """O(1) lookup from SVO filter ID string to row index."""
    return {str(f): i for i, f in enumerate(filter_ids)}
