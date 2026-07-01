#!/usr/bin/env python3
"""perf_report.py — P21 performance sweep on b5-vcatchment (Task 21.3.3).

Runs Frehg2 with synthetic scaling cases (cell count × step count), reads the per-run
simulation_summary.txt perf_* lines, and emits perf_report.{md,json} with a scaling table and
a simple roofline-style arithmetic-intensity estimate.

Usage:
    tools/perf_report.py --bin build/frehg2 --out build/perf
    tools/perf_report.py --quick   # small cases for CI
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# (label, nx, ny, max_steps) — native b5 is 101×5555 cells; synthetic flats for scaling.
FULL_CASES = [
    ("b5_native_100s", 101, 55, 100),
    ("cells_10k_1k", 100, 100, 1000),
    ("cells_100k_100s", 316, 316, 100),
    ("cells_1M_10s", 1000, 1000, 10),
]
QUICK_CASES = [
    ("quick_1k_50s", 32, 32, 50),
    ("quick_10k_20s", 100, 100, 20),
]


def parse_summary(path: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("frehg2") or line.startswith("simulation_id"):
                continue
            if line.startswith("modules"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            key, val = parts[0], parts[1].strip().strip('"')
            try:
                out[key] = float(val)
            except ValueError:
                pass
    return out


def make_yaml(nx: int, ny: int, max_steps: int, out_h5: str) -> str:
    import yaml

    doc = {
        "schema_version": "2.0",
        "simulation": {"id": "perf_sweep", "mode": "surface_water", "title": "P21 perf sweep"},
        "domain": {
            "nx": nx,
            "ny": ny,
            "nz": 1,
            "dx": 1.0,
            "dy": 1.0,
            "dz": 0.1,
            "botz": 0.0,
            "bathymetry": {"from_file": False, "constant": -0.05},
        },
        "time": {"dt": 1.0, "t_end": 1.0e9, "max_steps": max_steps, "output_interval": 0},
        "modules": {"surface_water": True, "groundwater": False, "solute": False},
        "surface_water": {
            "gravity": 9.81,
            "manning": 0.03,
            "min_depth": 1.0e-6,
            "viscosity": {"x": 0.0, "y": 0.0},
            "h_diffusion_ref": 0.1,
            "waterfall_depth": 1.0e-6,
        },
        "initial_conditions": {"surface_water": {"eta": 0.0}},
        "sources": {"surface": {"rainfall": {"from_file": False, "rate": 2.0e-5}, "evaporation": {"rate": 0.0}}},
        "output": {"format": "hdf5", "filename": out_h5, "io_mode": "serial_gather"},
    }
    return yaml.safe_dump(doc)


def run_case(binary: str, label: str, nx: int, ny: int, steps: int, out_dir: str) -> Dict[str, Any]:
    case_dir = os.path.join(out_dir, label)
    os.makedirs(case_dir, exist_ok=True)
    h5_rel = "out.h5"
    summary = os.path.join(case_dir, "simulation_summary.txt")
    try:
        import yaml  # noqa: F401
    except ImportError:
        return {"label": label, "error": "PyYAML required", "rc": 2}

    cfg_text = make_yaml(nx, ny, steps, h5_rel)
    cfg_path = os.path.join(case_dir, "run.yaml")
    with open(cfg_path, "w") as f:
        f.write(cfg_text)
    t0 = time.time()
    rc = subprocess.run([binary, os.path.abspath(cfg_path)], cwd=case_dir,
                        capture_output=True, text=True).returncode
    wall = time.time() - t0
    stats = parse_summary(summary)
    if rc != 0 and not stats:
        # Surface subprocess failure for CI/debugging.
        return {"label": label, "nx": nx, "ny": ny, "cells": nx * ny, "max_steps": steps,
                "rc": rc, "error": "frehg2 run failed", "wall_seconds": time.time() - t0}
    cells = nx * ny
    steps_done = stats.get("time_steps_completed", 0) or steps
    total_runtime = stats.get("total_runtime_seconds", wall)
    sw_ksp = stats.get("perf_sw_ksp_seconds", 0.0)
    sw_asm = stats.get("perf_sw_assembly_seconds", 0.0)
    sw_upd = stats.get("perf_sw_update_seconds", 0.0)
    ksp_iters = stats.get("perf_ksp_iterations", 0.0)
    bytes_staged = stats.get("perf_bytes_staged", 0.0)
    sec_per_step = total_runtime / max(steps_done, 1.0)
    # Roofline proxy: staged bytes / KSP time (GB/s effective staging bandwidth).
    bw = (bytes_staged / max(sw_ksp + sw_asm, 1e-12)) / 1.0e9
    return {
        "label": label,
        "nx": nx,
        "ny": ny,
        "cells": cells,
        "max_steps": steps,
        "steps_done": steps_done,
        "rc": rc,
        "wall_seconds": total_runtime,
        "seconds_per_step": sec_per_step,
        "perf_sw_ksp_seconds": sw_ksp,
        "perf_sw_assembly_seconds": sw_asm,
        "perf_sw_update_seconds": sw_upd,
        "perf_ksp_iterations": ksp_iters,
        "perf_bytes_staged": bytes_staged,
        "effective_staging_gbs": bw,
        "summary_path": summary,
    }


def scaling_quality(results: List[Dict[str, Any]]) -> List[str]:
    notes = []
    by_cells = sorted([r for r in results if r.get("rc") == 0], key=lambda r: r["cells"])
    for i in range(1, len(by_cells)):
        a, b = by_cells[i - 1], by_cells[i]
        cell_ratio = b["cells"] / max(a["cells"], 1)
        time_ratio = b["wall_seconds"] / max(a["wall_seconds"], 1e-12)
        linearity = time_ratio / cell_ratio
        notes.append(f"{a['label']}→{b['label']}: cell×{cell_ratio:.1f}, time×{time_ratio:.2f}, "
                     f"linearity={linearity:.2f} (1.0 = perfect)")
    return notes


def write_reports(out_dir: str, results: List[Dict[str, Any]], quick: bool) -> None:
    os.makedirs(out_dir, exist_ok=True)
    payload = {"quick": quick, "cases": results, "scaling_notes": scaling_quality(results)}
    with open(os.path.join(out_dir, "perf_report.json"), "w") as f:
        json.dump(payload, f, indent=2)

    lines = [
        "# Frehg2 Performance Report (P21)",
        "",
        f"Mode: {'quick smoke' if quick else 'full sweep'}",
        "",
        "| Case | cells | steps | wall [s] | s/step | sw_ksp [s] | sw_asm [s] | KSP iters | staging GB/s |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        if r.get("rc") != 0:
            lines.append(f"| {r['label']} | — | — | FAIL rc={r['rc']} | | | | | |")
            continue
        lines.append(
            f"| {r['label']} | {r['cells']} | {int(r['steps_done'])} | {r['wall_seconds']:.3f} | "
            f"{r['seconds_per_step']:.4e} | {r['perf_sw_ksp_seconds']:.3f} | "
            f"{r['perf_sw_assembly_seconds']:.3f} | {int(r['perf_ksp_iterations'])} | "
            f"{r['effective_staging_gbs']:.2f} |"
        )
    lines += ["", "## Scaling notes", ""]
    for n in payload["scaling_notes"]:
        lines.append(f"- {n}")
    lines += [
        "",
        "## Deferred gates (require Linux HPC / CUDA)",
        "",
        "- MPI parallel efficiency ≥ 0.7 at 16 ranks (Task 21.3.5) — not measured on macOS.",
        "- GPU speedup ≥ 5× vs OpenMP at 1M cells (Task 21.3.6) — deferred per gpu_validation_policy.md.",
        "",
    ]
    with open(os.path.join(out_dir, "perf_report.md"), "w") as f:
        f.write("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=os.path.join(REPO, "build", "frehg2"))
    ap.add_argument("--out", default=os.path.join(REPO, "build", "perf"))
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    binary = args.bin if os.path.isabs(args.bin) else os.path.join(REPO, args.bin)
    cases = QUICK_CASES if args.quick else FULL_CASES
    results = []
    for label, nx, ny, steps in cases:
        print(f"[perf] {label}: {nx}x{ny} × {steps} steps ...")
        results.append(run_case(binary, label, nx, ny, steps, args.out))
    write_reports(args.out, results, args.quick)
    failed = [r for r in results if r.get("rc") != 0]
    print(f"[perf] report -> {args.out}/perf_report.{{md,json}} ({len(failed)} failures)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
