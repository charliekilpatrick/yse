"""Reusable helpers for photometry I/O, passbands, plotting, HTTP, selection, bolometric export, and spectra."""

from yse.util.bands import uniformize_yse_xy_passbands
from yse.util.bolometric import (
    EXTRABOL_FLT_TO_SVO,
    SUPERBOL_FLT_REPLACEMENTS,
    extrabol_format_table,
    snana_passband_extrabol_table,
    write_extrabol_dat,
    write_superbol_filter_files,
)
from yse.util.http import (
    ensure_dir,
    env_optional,
    env_required,
    fetch_bytes,
    fetch_text,
    install_basic_auth,
    verify_size_or_remove,
)
from yse.util.io import (
    list_snana_txt,
    load_thesis_yse_photometry_tables,
    read_whitespace_lightcurve,
    read_yse_ascii_photometry,
    unique_filters,
)
from yse.util.logging_config import configure_logging
from yse.util.plotting import (
    THESIS_BAND_COLORS,
    plot_dr1_snana_bands,
    plot_multiband_lightcurve,
    plot_normalized_spectra_csv,
)
from yse.util.selection import (
    griz_pre_post_peak_counts,
    peak_mjd_from_min_mag,
    select_snana_by_griz_depth,
)
from yse.util.spectrum import smooth_flux_300kms

__all__ = [
    "configure_logging",
    "EXTRABOL_FLT_TO_SVO",
    "SUPERBOL_FLT_REPLACEMENTS",
    "THESIS_BAND_COLORS",
    "ensure_dir",
    "env_optional",
    "env_required",
    "extrabol_format_table",
    "fetch_bytes",
    "griz_pre_post_peak_counts",
    "fetch_text",
    "install_basic_auth",
    "list_snana_txt",
    "load_thesis_yse_photometry_tables",
    "peak_mjd_from_min_mag",
    "read_whitespace_lightcurve",
    "read_yse_ascii_photometry",
    "plot_dr1_snana_bands",
    "plot_multiband_lightcurve",
    "plot_normalized_spectra_csv",
    "select_snana_by_griz_depth",
    "smooth_flux_300kms",
    "snana_passband_extrabol_table",
    "uniformize_yse_xy_passbands",
    "unique_filters",
    "verify_size_or_remove",
    "write_extrabol_dat",
    "write_superbol_filter_files",
]
