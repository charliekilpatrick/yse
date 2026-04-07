# Command-line tools (`yse.cli`)

```bash
python -m yse.cli.<module> [args...]
```

Shared code: **`yse.util`** (photometry I/O, plotting, passbands, bolometric export, HTTP, selection, spectra).

## Analysis

| Module | Description |
|--------|-------------|
| **`snana_cuts`** | DR1 SNANA Type II listing, cuts, metadata. |
| **`snana_depth`** | DR1 SNANA band plots + YSE-PZ $griz$ depth selection. |
| **`thesis_phot`** | Thesis `all_YSE_phot`: tables, plots, extrabol / superbol inputs. |
| **`nickel_overview`** | 2020jfo multiband photometry + bolometric comparison. |
| **`nickel_mcmc`** | Valenti-style bolometric fit + emcee + $^{56}$Ni vs SN 1987A. |
| **`host_plots`** | Literature + BLAST SN–host correlation figures. Default output: `yse.paths.NICKEL_ANALYSIS_FIGURES`. |
| **`superbol_plots`** | Batch bolometric / multiband plots from superbol outputs. |
| **`smooth_spectrum`** | Spectrum smoothing helper. |
| **`extrabol_helpers`** | Extrabol `.dat` filter checks; bolometric ascii helpers. |
| **`spectrum_plots`** | Normalized spectra from CSVs in a working directory. |
| **`villar_gp`** | II-P Villar/Bazin GP fit per band (needs **light-curve**, in base deps). |
| **`ghost_demo`** | Minimal GHOST example — `pip install -e ".[ghost]"`. |

## Downloads and maintenance

| Module | Description |
|--------|-------------|
| **`lick_download`** | UCOLICK / Lick APF, Kast, or Nickel FITS fetch (`yse/config/env.example`). |
| **`vlass_cutouts`** | VLASS quicklook image lookup (demo coordinates). |
| **`flatten_data`** | Migrates old `data/yse_dr1/...` tree into flat `yse/data/`. |
| **`jupytext_export`** | Extracts code cells from `.ipynb` into `yse/notebook_cell_exports/`. |

## GP extrabol (`yse.extrabol`)

```bash
python -m yse.extrabol.extrabol_gp --help
```

Nebular SED features require **`pip install -e ".[extrabol-gp]"`** (synphot + numba).
