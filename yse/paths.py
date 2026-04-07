"""Repository root and locations under ``yse/data/``."""

from __future__ import annotations

from pathlib import Path

# Repository root (parent of the ``yse`` package directory)
REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "yse" / "data"
DR1_SNANA_II = DATA / "dr1_snana_ii"
DR1_SNANA_NOT_II = DATA / "dr1_snana_not_ii"
THESIS = DATA / "thesis"
THESIS_PHOT = THESIS / "all_YSE_phot"
NICKEL_DIR = THESIS_PHOT / "nickel"
# Valenti / Anderson / BLAST JSON / YSE nickel tables (thesis tree)
NICKEL_DATA = NICKEL_DIR / "data"
# Default PNG/CSV output for ``yse.cli.host_plots`` (created on demand)
NICKEL_ANALYSIS_FIGURES = NICKEL_DIR / "figures" / "nickel_analysis"
SUPERBOL_ROOT = THESIS_PHOT / "superbol-master"
EXTRABOL_INPUTS = THESIS_PHOT / "extrabol_inputs"
YSEPZ = DATA / "ysepz"
IIB_FLASH = DATA / "iib_iin_flash"
IIP_VILLAR = DATA / "iip_villar"
SPECTRA = DATA / "yse_spectra"
# GP extrabol pipeline (synphot + optional Numba)
EXTRABOL = REPO / "yse" / "extrabol"
