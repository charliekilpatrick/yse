#!/usr/bin/env python3
"""One-shot migrator: old ``data/yse_dr1/...`` → flat ``data/`` layout."""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "yse" / "data"


def main() -> None:
    ydr = DATA / "yse_dr1"
    if not ydr.is_dir():
        print("No data/yse_dr1 to flatten.", file=sys.stderr)
        return
    type_ii = ydr / "type_ii"
    moves: list[tuple[Path, Path]] = [
        (type_ii / "snana", DATA / "dr1_snana_ii"),
        (ydr / "not_type_ii", DATA / "dr1_snana_not_ii"),
        (type_ii / "ysepz", DATA / "ysepz"),
        (type_ii / "iib_iin_flash", DATA / "iib_iin_flash"),
        (type_ii / "thesis_sample", DATA / "thesis"),
    ]
    for src, dst in moves:
        if not src.exists():
            print(f"Skip missing {src}", file=sys.stderr)
            continue
        if dst.exists():
            print(f"Refuse overwrite {dst}", file=sys.stderr)
            sys.exit(1)
        shutil.move(str(src), str(dst))
    shutil.rmtree(ydr, ignore_errors=True)
    iv = DATA / "iip_villar_fit"
    if iv.is_dir() and not (DATA / "iip_villar").exists():
        shutil.move(str(iv), str(DATA / "iip_villar"))
    print("Flatten complete.")


if __name__ == "__main__":
    main()
