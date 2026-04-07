"""VLASS quicklook image lookup and download (NRAO archive)."""

from __future__ import annotations

import re
import shutil

import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.utils.data import download_file


def get_tile_definitions() -> Table:
    url = "https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php"
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    data = [d.split() for d in r.text.split("\n") if d]
    data = data[3:]
    data = list(map(list, zip(*data)))
    names = ("tile", "dmin", "dmax", "rmin", "rmax", "available")
    tab = Table(data, names=names)
    for col in ("dmin", "dmax", "rmin", "rmax"):
        tab[col] = tab[col].astype(float)
    return tab


def get_tile(ra: float, dec: float, tab: Table):
    for row in tab:
        if ra > row["rmin"] and ra < row["rmax"] and dec > row["dmin"] and dec < row["dmax"]:
            return row["tile"]
    return None


def download_image(ra, dec, obj, tile, tab):
    if not tile:
        return None

    row = tab[tab["tile"] == tile][0]
    objcoord = SkyCoord(ra, dec, unit="deg")
    base_url = "https://archive-new.nrao.edu/vlass/quicklook/{0}/{1}/"
    base_url = base_url.format(row["available"], tile)
    r = requests.get(base_url, timeout=120)
    r.raise_for_status()
    images = [d for d in r.text.split("\n") if ("<a href" in d and "VLASS" in d)]
    img_data = []
    for im in images:
        im = im.split(">")[1]
        im = im.split("<")[0]
        if im.strip():
            coord_data = im.split(".")[4]
            r_hms = coord_data[1:7]
            d_dms = coord_data[7:15]
            ra_hms = r_hms[0:2] + ":" + r_hms[2:4] + ":" + r_hms[4:6]
            de_dms = d_dms[0:3] + ":" + d_dms[3:5] + ":" + d_dms[5:7]
            coord = SkyCoord(ra_hms, de_dms, unit=(u.hour, u.deg))
            img_data.append((im, coord))

    img_data = sorted(img_data, key=lambda x: objcoord.separation(x[1]).degree)
    closest = img_data[0]
    url = base_url + closest[0]
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    s = re.search(r"VLASS.*?\.tt0\.subim\.fits", r.text)
    if not s:
        return None
    url += "/" + s[0]
    dat = download_file(url)
    shutil.move(dat, obj + ".fits")
    return url


def main() -> None:
    tab = get_tile_definitions()
    data = [["23:37:16.46", "53:24:57.7"]]
    for ra, dec in data:
        coord = SkyCoord(ra, dec, unit=(u.hour, u.deg))
        rh = coord.ra.hour
        rd = coord.ra.degree
        tile = get_tile(rh, coord.dec.degree, tab)
        download_image(rd, coord.dec.degree, str(round(coord.dec.degree)), tile, tab)


if __name__ == "__main__":
    main()
