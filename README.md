# YSE Type II population

YSE DR1 Type II supernovae: SNANA light curves, thesis photometry, and Python tooling.

## Where things live

- **`yse/data/`** — Datasets only (see `yse/data/README.md`).
- **`yse/paths.py`**, **`yse/snana_io.py`** — Repository data locations and SNANA light-curve reader (re-exported from `import yse`).
- **`yse/population/`** — Project pipelines: **`analysis/`** (CLIs), **`downloads/`** (UCOLICK / Lick + VLASS fetch), **`notebooks/`** (a few consolidated exploratory CLIs; see `notebooks/README.md`), **`utils/`** (shared helpers including HTTP used by downloads), **`devtools/`** (migrators / nb export).
- **`yse/population/analysis/`** — Analysis CLIs plus bundled *extrabol* / *superbol* legacy scripts (`extrabol.py`, `superbol_interactive.py`) and `vendor_tools.py` (spectrum smoothing).
- **`yse/config/env.example`** — Environment variable names for secrets and machine paths (copy to a private file; never commit real values).

There are no Jupyter notebooks in the repo; no plaintext passwords or usernames should be committed.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -U pip setuptools wheel
pip install -e .
pip install -e ".[thesis]"    # emcee — nickel fits
```

## Quick analysis

```bash
python -m yse.population.analysis.snana_population cuts --print-names
python -m yse.population.analysis.nickel_valenti_emcee rth
```

See `yse/docs/README.md` for download / cron scripts (require `yse/config/env.example` variables).

## License

MIT unless survey policies restrict redistribution of data files.
