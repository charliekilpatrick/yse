# yse-type-ii-population

Python package for **Young Supernova Experiment DR1** spectroscopic Type II SNe: SNANA light-curve I/O, thesis photometry pipelines, bolometric and \(^{56}\)Ni analysis, BLAST host–SN figures, and download helpers.

**Python 3.10+** required. Optional nebular SED code in `yse.extrabol.extrabol_gp` needs **`pip install -e ".[extrabol-gp]"`** (synphot + numba).

## Layout

| Path | Role |
|------|------|
| **`yse/paths.py`** | Canonical paths under `yse/data/`. |
| **`yse/snana.py`** | SNANA `*.dat` directory reader (also `from yse import read_YSE_ZTF_snana_dir`). |
| **`yse/util/`** | Shared I/O, plotting, passbands, bolometric export, HTTP, selection, spectra. |
| **`yse/cli/`** | Commands: `python -m yse.cli.<name>`. See **`yse/cli/README.md`**. |
| **`yse/extrabol/`** | Gaussian-process extrabol (`extrabol_gp`). |
| **`yse/data/`** | Data only — not importable (see **`yse/data/README.md`**). |
| **`yse/config/`** | Environment templates (no secrets). |
| **`yse/docs/`** | Notes and background reading. |

No notebooks are tracked; optional exports go to `yse/notebook_cell_exports/` via `yse.cli.jupytext_export`.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -U pip setuptools wheel
pip install -e .
```

```bash
pip install -e ".[extrabol-gp]"   # synphot + numba — nebular extrabol_gp
pip install -e ".[ghost]"          # astro-ghost
pip install -e ".[spectroscopy]"   # pypeit
pip install -e ".[google]"         # Google API clients
pip install -e ".[dev]"            # pytest + ruff
```

## Quick commands

```bash
python -m yse.cli.snana_cuts cuts --print-names
python -m yse.cli.nickel_mcmc rth
python -m yse.cli.host_plots
```

Default PNG output for `host_plots`: `yse.paths.NICKEL_ANALYSIS_FIGURES`.

## Development

```bash
pip install -e ".[dev]"
ruff check yse
```

## License

MIT unless archive policies restrict redistribution of files under `yse/data/`.

## Citation

Cite YSE DR1 and surveys you use; acknowledge this software as appropriate for your journal.
