"""HTTP helpers: env lookups, basic auth, and small urllib fetches (used by archive download CLIs)."""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path


def env_required(key: str) -> str:
    v = os.environ.get(key, "").strip()
    if not v:
        raise RuntimeError(
            f"Missing required environment variable {key}. See yse/config/env.example."
        )
    return v


def env_optional(key: str, default: str | None = None) -> str | None:
    v = os.environ.get(key, "").strip()
    return v if v else default


def install_basic_auth(base_url: str, username: str, password: str) -> None:
    passman = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    passman.add_password(None, base_url, username, password)
    auth = urllib.request.HTTPBasicAuthHandler(passman)
    opener = urllib.request.build_opener(auth)
    urllib.request.install_opener(opener)


def fetch_bytes(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def fetch_text(url: str, timeout: int = 120) -> str:
    return fetch_bytes(url, timeout=timeout).decode("utf-8", errors="replace")


def verify_size_or_remove(path: Path, expected: int) -> bool:
    if expected <= 0:
        return True
    got = path.stat().st_size
    if got != expected:
        print(f"WARNING: {path.name} has wrong size ({got}); deleting.")
        path.unlink(missing_ok=True)
        return False
    return True


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)
