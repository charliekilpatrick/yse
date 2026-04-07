# Data (flat layout)

All inputs and generated tables stay **out of** the Python package tree. Top-level buckets under `yse/data/`:

| Directory | Role |
|-----------|------|
| `dr1_snana_ii/` | YSE DR1 Type II SNANA `*.dat` light curves. |
| `dr1_snana_not_ii/` | DR1 objects not broad Type II. |
| `thesis/` | Thesis sample under `all_YSE_phot/`: YSE-PZ photometry (`*.snana.txt`), per-object SNANA-style `2020*.dat` dumps, `extrabol_inputs/`, superbol outputs, nickel tables, etc. Use `yse.paths.THESIS_PHOT` (alias `THESIS_SNANA_DAT` for the flat `*.dat` location). |
| `ysepz/` | YSE-PZ exports for the Type II context. |
| `iib_iin_flash/` | IIb / IIn / flash subsample. |
| `extrabol_19mhm/`, `extrabol_inputs/` | Standalone extrabol test / input trees. |
| `iip_villar/` | II-P Villar-style fit inputs and outputs. |
| `yse_spectra/`, `yse_spectra.tar` | Reduced spectra. |

Deeper folders under `thesis/all_YSE_phot/` keep filenames from colliding. Code: `yse.util` (shared helpers), `yse.cli` (commands), `yse.extrabol` (GP extrabol), plus `yse.paths` and `yse.snana` (see the repository `README.md`).
