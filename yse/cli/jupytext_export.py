#!/usr/bin/env python3
"""Extract code cells from every .ipynb under the repo into ``yse/notebook_cell_exports/``."""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "yse" / "notebook_cell_exports"


def sanitize(name: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9_]+", "_", name).strip("_")
    return s or "notebook"


def export_one(ipynb: Path) -> Path:
    rel = ipynb.relative_to(ROOT)
    with ipynb.open(encoding="utf-8") as f:
        nb = json.load(f)
    chunks: list[str] = [
        f'"""Auto-exported from {rel.as_posix()} — edit in Python only; notebook removed."""\n',
        "# flake8: noqa\n",
    ]
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        src = cell.get("source", [])
        text = src if isinstance(src, str) else "".join(src)
        if not text.strip():
            continue
        chunks.append("\n")
        chunks.append(text)
        if not text.endswith("\n"):
            chunks.append("\n")
    stem = sanitize(ipynb.stem)
    parent = sanitize(ipynb.parent.name)
    out_name = f"{parent}__{stem}.py" if parent else f"{stem}.py"
    out_path = OUT / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(chunks), encoding="utf-8")
    return out_path


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    seen: dict[str, int] = {}
    for ipynb in sorted(ROOT.rglob("*.ipynb")):
        if ".git" in ipynb.parts:
            continue
        out = export_one(ipynb)
        key = out.name
        if key in seen:
            seen[key] += 1
            alt = out.with_name(f"{out.stem}_{seen[key]}{out.suffix}")
            out.rename(alt)
            out = alt
        else:
            seen[key] = 0
        print(out.relative_to(ROOT))


if __name__ == "__main__":
    main()
