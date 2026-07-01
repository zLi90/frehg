#!/usr/bin/env python3
"""Generate Frehg2 b3-kirkland per-layer soil-class files from the legacy SERGHEI soilID.input.

The SERGHEI `soilID.input` is a 2-D (z-rows x x-cols) class grid for the 50 x 30 (x x z)
Kirkland domain (ny = 1). Frehg2's 3-D soil map (`soil.map.layers`) wants one class file per
vertical layer (top k=0 -> bottom k=nz-1), each holding nx*ny class ids in row-major
gi + gj*nx order. soilID row r (top-first) maps directly to Frehg2 layer k=r.

Run from the repo root (or anywhere): writes benchmarks/b3-kirkland/soil_layers/layer_kk.txt.
"""
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent.parent / "legacy" / "benchmarks" / "b3-kirkland" / "input" / "soilID.input"
OUT = HERE / "soil_layers"


def main() -> None:
    lines = SRC.read_text().strip().splitlines()
    # First line is "n_soil 2"; the rest are the 30 z-rows of 50 class ids.
    rows = [ln.split() for ln in lines[1:] if ln.strip()]
    nz = len(rows)
    nx = len(rows[0])
    assert nz == 30 and nx == 50, f"unexpected soilID shape {nz}x{nx}"
    OUT.mkdir(exist_ok=True)
    for k, row in enumerate(rows):
        assert len(row) == nx, f"row {k} has {len(row)} cols"
        (OUT / f"layer_{k:02d}.txt").write_text(" ".join(row) + "\n")
    print(f"wrote {nz} layer files (nx={nx}, ny=1) to {OUT}")


if __name__ == "__main__":
    main()
