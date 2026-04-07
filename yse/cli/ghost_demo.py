#!/usr/bin/env python3
"""Minimal GHOST host association example (requires ``pip install astro-ghost``)."""
# flake8: noqa

from __future__ import annotations

from pathlib import Path
from astropy.coordinates import SkyCoord
from astro_ghost.ghostHelperFunctions import getGHOST, getTransientHosts

# real=False avoids downloading the full transient DB; set real=True for production.
getGHOST(real=False, verbose=1)

transient_position = SkyCoord(ra=160.9708, dec=11.6714, unit="deg")
host_data = getTransientHosts(
    transientCoord=[transient_position],
    transientName=["2012aw"],
    verbose=1,
    starcut="gentle",
    ascentMatch=False,
)
print(host_data)

for sub in ("hostSpectra", "SNspectra", "hostPostageStamps"):
    p = Path(f"./{sub}/")
    p.mkdir(parents=True, exist_ok=True)
