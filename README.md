# YSE Type II population

Python package and CLIs for the YSE DR1 Type II supernova sample: SNANA light curves, thesis photometry, nickel mass work, and supporting downloads. **Requires Python 3.10+** (needed for optional [synphot](https://synphot.readthedocs.io/)–based nebular code in `extrabol_gp`).

## Layout

| Location | Role |
|----------|------|
| **`yse/paths.py`** | Canonical `Path` objects under `yse/data/` (thesis tree, DR1 SNANA dirs, nickel, etc.). |
| **`yse/snana_io.py`** | Reader for YSE/ZTF SNANA `*.dat` directories (re-exported from `import yse`). |
| **`yse/lib/`** | Shared library code: photometry I/O, passbands, plotting, selection, HTTP helpers, bolometric export. |
| **`yse/cli/`** | Runnable modules (`python -m yse.cli.<name>`). See **`yse/cli/README.md`** for the full command list. |
| **`yse/extrabol/`** | GP **`extrabol_gp`** pipeline ([synphot](https://synphot.readthedocs.io/) + optional Numba). |
| **`yse/data/`** | Datasets only — not importable (see **`yse/data/README.md`**). |
| **`yse/config/`** | Example env / account templates (no secrets committed). |
| **`yse/docs/`** | Non-code notes and background reading. |

There are no Jupyter notebooks in the repository. Auto-exports from stray `.ipynb` files (if you run the exporter) land in **`yse/notebook_cell_exports/`**.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -U pip setuptools wheel
pip install -e .
```

Optional stacks (only if you need those tools):

```bash
pip install -e ".[extrabol-gp]"   # synphot + numba — nebular path in extrabol_gp
pip install -e ".[ghost]"         # astro-ghost — host cross-match demo
pip install -e ".[spectroscopy]"  # pypeit
pip install -e ".[google]"        # Google API clients
```

Core analysis dependencies (**emcee**, **george**, **astroquery**, **extinction**, **light-curve**) are required by the default install so every CLI in `yse.cli` resolves without extra extras, except where noted in **`yse/cli/README.md`**.

## Quick commands

```bash
python -m yse.cli.snana_type_ii_cuts cuts --print-names
python -m yse.cli.nickel_mcmc_valenti rth
python -m yse.cli.nickel_host_blast_figures
```

### Migrating from the old layout

The former `yse.population.utils` package is now **`yse.lib`**. Runnable modules moved from `yse.population.analysis`, `notebooks`, `downloads`, and `devtools` into **`yse.cli`** with clearer names (for example `snana_population` → `snana_type_ii_cuts`, `nickel_host_figures` → `nickel_host_blast_figures`). GP extrabol lives under **`yse.extrabol`**. Update any `python -m yse.population…` invocations using the table in **`yse/cli/README.md`**.

UCOLICK / Lick archive downloads use variables described in **`yse/config/env.example`**. See **`yse/docs/README.md`** for operational notes.

## License

MIT unless survey or archive policies restrict redistribution of bundled data files.
