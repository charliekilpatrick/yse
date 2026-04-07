#!/usr/bin/env python3
"""UCOLICK / Lick archive downloads: APF, Kast (Shane), or Nickel FITS."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from astropy.time import Time
from dateutil.parser import parse

from yse.util.http import (
    ensure_dir,
    env_optional,
    env_required,
    fetch_bytes,
    fetch_text,
    install_basic_auth,
    verify_size_or_remove,
)

UCOLICK_BASE = "https://mthamilton.ucolick.org"
APF_EXPECTED_SIZE = 19189440


def _build_apf_url(date_yyyy_mm_dd: str, archive_user_folder: str) -> str:
    return f"{UCOLICK_BASE}/data/{date_yyyy_mm_dd}/APF/{archive_user_folder}/"


def _build_kast_url(date_yyyy_mm_dd: str, archive_user_folder: str) -> str:
    return f"{UCOLICK_BASE}/data/{date_yyyy_mm_dd}/shane/{archive_user_folder}/"


def _kast_expected_size(name: str) -> int:
    if name.startswith("b"):
        return 1497600
    if name.startswith("r"):
        return 3968640
    return 0


def download_apf(
    download_date,
    rawdir: Path,
    *,
    http_user: str,
    http_password: str,
    archive_user_folder: str,
) -> int:
    date_s = download_date.strftime("%Y-%m-%d")
    base_url = _build_apf_url(date_s, archive_user_folder)
    install_basic_auth(base_url, http_user, http_password)
    try:
        html = fetch_text(base_url)
    except OSError as e:
        print(f"No data for {date_s}: {e}", file=sys.stderr)
        return 0
    dest_root = rawdir / download_date.strftime("ut%y%m%d")
    ensure_dir(dest_root)
    n = 0
    for name in re.findall(r'value="(.+\.fits)"', html):
        local = dest_root / name
        url = base_url + name
        if not local.is_file():
            print(f"Downloading {date_s}: {local}")
            local.write_bytes(fetch_bytes(url))
            n += 1
        if not verify_size_or_remove(local, APF_EXPECTED_SIZE):
            continue
        print(f"OK {local}")
    return n


def download_kast(
    download_date,
    rawdir: Path,
    *,
    http_user: str,
    http_password: str,
    archive_user_folder: str,
) -> None:
    date_s = download_date.strftime("%Y-%m-%d")
    base_url = _build_kast_url(date_s, archive_user_folder)
    install_basic_auth(base_url, http_user, http_password)
    try:
        html = fetch_text(base_url)
    except OSError as e:
        print(f"No Kast data for {date_s}: {e}", file=sys.stderr)
        sys.exit(1)
    dest_root = rawdir / download_date.strftime("ut%y%m%d")
    ensure_dir(dest_root)
    for name in re.findall(r'value="([br].+\.fits)"', html):
        local = dest_root / name
        url = base_url + name
        exp = _kast_expected_size(name)
        if not local.is_file():
            print(f"Downloading {local}")
            local.write_bytes(fetch_bytes(url))
            try:
                local.chmod(0o775)
            except OSError:
                pass
        if exp and not verify_size_or_remove(local, exp):
            continue
        print(f"OK {local}")


def load_nickel_accounts(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, list):
        raise ValueError("Nickel accounts JSON must be a list")
    out = []
    for row in data:
        if not isinstance(row, dict):
            continue
        af, hu, hp = row.get("archive_folder"), row.get("http_user"), row.get("http_password")
        if af is None or hu is None or hp is None:
            continue
        out.append(
            {"archive_folder": str(af), "http_user": str(hu), "http_password": str(hp)}
        )
    if not out:
        raise ValueError("No valid entries in accounts JSON")
    return out


def _try_nickel_page(nickel_date: str, account: dict) -> tuple[str, str] | None:
    folder = account["archive_folder"]
    nickel_url = f"{UCOLICK_BASE}/data/{nickel_date}/nickel/{folder}"
    install_basic_auth(nickel_url, account["http_user"], account["http_password"])
    try:
        return nickel_url, fetch_text(nickel_url, timeout=30)
    except OSError:
        return None


def download_nickel(stagedir: Path, date_arg: str, accounts_path: Path) -> None:
    nickel_date = Time(date_arg).datetime.strftime("%Y-%m-%d")
    accounts = load_nickel_accounts(accounts_path)
    result = None
    for acct in accounts:
        result = _try_nickel_page(nickel_date, acct)
        if result:
            break
    if result is None:
        print("Could not open any Nickel archive URL for this date.", file=sys.stderr)
        sys.exit(1)
    nickel_url, html = result
    print(f"Using {nickel_url}")
    ensure_dir(stagedir)
    names = [f"d{n}.fits" for n in re.findall(r'value="d(.+)\.fits"', html)]
    for name in names:
        full = stagedir / name
        if full.is_file():
            print(f"{full} already exists")
            continue
        url = f"{nickel_url}/{name}"
        print(f"Downloading {url}")
        full.write_bytes(fetch_bytes(url))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    pa = sub.add_parser("apf", help="Download APF FITS for one night")
    pa.add_argument("date", help="Observation date (dateutil)")
    pa.add_argument("rawdir", nargs="?", type=Path, help="Base raw dir or env UCOLICK_RAWDIR_APF")
    pa.set_defaults(func=_cli_apf)

    pk = sub.add_parser("kast", help="Download Kast FITS for one night")
    pk.add_argument("date", help="Observation date")
    pk.add_argument("rawdir", nargs="?", type=Path, help="Base raw dir or env UCOLICK_RAWDIR_KAST")
    pk.set_defaults(func=_cli_kast)

    pn = sub.add_parser("nickel", help="Download Nickel FITS for one night")
    pn.add_argument("stagedir", type=Path, help="Output directory")
    pn.add_argument("date", help="Date string for astropy Time")
    pn.add_argument("--accounts", type=Path, help="Accounts JSON (default env)")
    pn.set_defaults(func=_cli_nickel)

    args = p.parse_args()
    args.func(args)


def _cli_apf(args: argparse.Namespace) -> None:
    dt = parse(args.date)
    raw = args.rawdir or Path(env_optional("UCOLICK_RAWDIR_APF") or Path.cwd() / "apf_raw")
    download_apf(
        dt,
        raw.expanduser().resolve(),
        http_user=env_required("UCOLICK_HTTP_USER"),
        http_password=env_required("UCOLICK_HTTP_PASSWORD"),
        archive_user_folder=env_required("UCOLICK_APF_ARCHIVE_USER"),
    )


def _cli_kast(args: argparse.Namespace) -> None:
    dt = parse(args.date)
    raw = args.rawdir or Path(env_optional("UCOLICK_RAWDIR_KAST") or Path.cwd() / "kast_raw")
    download_kast(
        dt,
        raw.expanduser().resolve(),
        http_user=env_required("UCOLICK_HTTP_USER"),
        http_password=env_required("UCOLICK_HTTP_PASSWORD"),
        archive_user_folder=env_required("UCOLICK_KAST_ARCHIVE_USER"),
    )


def _cli_nickel(args: argparse.Namespace) -> None:
    acc = args.accounts or Path(env_required("UCOLICK_NICKEL_ACCOUNTS_JSON"))
    download_nickel(
        args.stagedir.expanduser().resolve(),
        args.date,
        acc.expanduser().resolve(),
    )


if __name__ == "__main__":
    main()
