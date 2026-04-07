"""GP extrabol pipeline code kept separate from :mod:`yse.cli`.

``extrabol_gp`` uses :mod:`yse.extrabol._extrabol_accel` for faster blackbody / chi^2 paths,
and :mod:`yse.extrabol._synphot_phot` for nebular photometry via `synphot` (not pysynphot).
Optional ``numba`` is listed under the ``extrabol-gp`` extra in ``pyproject.toml``.
"""
