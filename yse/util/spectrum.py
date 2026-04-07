"""1D spectral smoothing (legacy 300 km/s Gaussian kernel on binned data)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve


def smooth_flux_300kms(
    path_in: str | Path,
    path_out: str | Path | None = None,
    *,
    sigma_bins: float = 1.28,
) -> Path:
    """
    Convolve column 1 of a two-column ASCII spectrum with a Gaussian kernel.

    Default ``sigma_bins=1.28`` matches 300 km/s FWHM when bins are 100 km/s (FWHM ≈ 2.35σ).
    """
    path_in = Path(path_in)
    if path_out is None:
        path_out = path_in.with_name(path_in.stem + "_smoothed300kms.txt")
    else:
        path_out = Path(path_out)

    data = np.loadtxt(path_in)
    kernel = Gaussian1DKernel(sigma_bins)
    y_smooth = convolve(data[:, 1], kernel)
    out = np.c_[data[:, 0], y_smooth]
    np.savetxt(path_out, out, fmt="%12.4e")
    return path_out
