"""
Synthetic photometry via `synphot` (https://synphot.readthedocs.io/en/latest/ ),
replacing legacy ``pysynphot`` in :mod:`yse.extrabol.extrabol_gp`.
"""

from __future__ import annotations

import numpy as np
import astropy.units as u
from synphot import Observation, SourceSpectrum, SpectralElement
from synphot.models import Empirical1D
from synphot import units as su
from synphot.spectrum import lazy_load_vega

_vega_cache: object | None = None


def vega_reference() -> object:
    global _vega_cache
    if _vega_cache is None:
        _vega_cache = lazy_load_vega()
    return _vega_cache


def source_spectrum_flam(wave_angstrom: np.ndarray, flux_flam: np.ndarray) -> SourceSpectrum:
    """Empirical source spectrum; ``wave`` in Å, ``flux`` in FLAM (erg/s/cm²/Å)."""
    wave = np.asarray(wave_angstrom, dtype=np.float64)
    flux = np.asarray(flux_flam, dtype=np.float64)
    return SourceSpectrum(
        Empirical1D,
        points=wave * u.AA,
        lookup_table=flux * su.FLAM,
    )


def spectral_element_throughput(
    wave_angstrom: np.ndarray,
    throughput: np.ndarray,
    *,
    name: str | None = None,
) -> SpectralElement:
    """Bandpass from wavelength (Å) and dimensionless throughput (0–1)."""
    wave = np.asarray(wave_angstrom, dtype=np.float64)
    thr = np.asarray(throughput, dtype=np.float64)
    kw: dict = dict(
        points=wave * u.AA,
        lookup_table=thr * u.dimensionless_unscaled,
    )
    if name is not None:
        kw["name"] = name
    return SpectralElement(Empirical1D, **kw)


def observation_effective_magnitude(
    sp: SourceSpectrum,
    bp: SpectralElement,
    binset_angstrom: np.ndarray,
    phot_unit: str,
) -> float:
    """
    Effective magnitude through ``bp`` (matches pysynphot ``Observation.effstim``).

    Parameters
    ----------
    phot_unit
        ``'AB'`` or ``'vega'`` (case-insensitive).
    """
    binset = np.asarray(binset_angstrom, dtype=np.float64)
    obs = Observation(sp, bp, binset=binset)
    pu = phot_unit.lower().strip()
    if pu == "ab":
        return float(obs.effstim(u.ABmag).value)
    if pu in ("vega", "vegamag"):
        return float(obs.effstim(su.VEGAMAG, vegaspec=vega_reference()).value)
    raise ValueError(f"Unsupported photometry unit {phot_unit!r}; use 'AB' or 'vega'.")
