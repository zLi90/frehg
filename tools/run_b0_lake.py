#!/usr/bin/env python3
"""run_b0_lake.py — well-balanced "lake at rest" check (P19 Task 19.3.6).

b0-lake is a SEPARATE check, NOT part of the b1-b6 suite: a flat free surface (eta = 1.0) over a
Gaussian bump must stay at rest. PASS criterion: max|eta - eta0| <= 1e-12 (machine precision)
across all output snapshots. (Velocity preservation is additionally asserted by the C++ ctest
test_swe_b0; here we check the eta well-balanced invariant end-to-end through the driver + HDF5.)

Usage: run_b0_lake.py [--bin build/frehg2] [--out <dir>] [--tol 1e-12]
Returns a JSON-able dict via run_b0(); exits non-zero on FAIL when run as a script.
"""
import argparse
import math
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ETA0 = 1.0


def _abs(p):
    return p if os.path.isabs(p) else os.path.join(REPO, p)


def run_b0(binary="build/frehg2", out_dir=None, tol=1e-12, max_steps=None):
    import h5py  # local import so importing this module never hard-requires h5py
    cfg = os.path.join(REPO, "benchmarks", "b0-lake", "b0-lake.yaml")
    if not os.path.exists(cfg):
        return {"id": "b0-lake", "tier": "FAIL", "reason": "config missing", "max_eta_dev": None}

    run_cfg = cfg
    tmp = None
    out_dir = out_dir or tempfile.mkdtemp(prefix="b0lake_")
    os.makedirs(out_dir, exist_ok=True)
    h5_path = os.path.join(out_dir, "b0-lake.h5")
    try:
        import yaml
        with open(cfg) as f:
            doc = yaml.safe_load(f)
        doc.setdefault("output", {})["filename"] = h5_path
        if max_steps is not None:
            doc.setdefault("time", {})["max_steps"] = int(max_steps)
        tmp = tempfile.NamedTemporaryFile("w", suffix=".yaml",
                                          dir=os.path.dirname(cfg), delete=False)
        yaml.safe_dump(doc, tmp)
        tmp.close()
        run_cfg = tmp.name

        rc = subprocess.run([_abs(binary), run_cfg], cwd=os.path.dirname(cfg)).returncode
        if rc != 0:
            return {"id": "b0-lake", "tier": "FAIL", "reason": f"run rc={rc}", "max_eta_dev": None}

        max_dev = 0.0
        finite = True
        with h5py.File(h5_path, "r") as f:
            grp = f.get("/surface/eta")
            if grp is None:
                return {"id": "b0-lake", "tier": "FAIL", "reason": "no /surface/eta",
                        "max_eta_dev": None}
            for key in grp.keys():
                arr = grp[key][:]
                for v in arr:
                    if not math.isfinite(float(v)):
                        finite = False
                    else:
                        max_dev = max(max_dev, abs(float(v) - ETA0))
        tier = "PASS" if (finite and max_dev <= tol) else "FAIL"
        return {"id": "b0-lake", "tier": tier, "max_eta_dev": max_dev, "tol": tol,
                "finite": finite,
                "reason": "" if tier == "PASS" else f"max|eta-1|={max_dev:.3e} > {tol:.0e}"}
    finally:
        if tmp and os.path.exists(tmp.name):
            os.unlink(tmp.name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="build/frehg2")
    ap.add_argument("--out", default=None)
    ap.add_argument("--tol", type=float, default=1e-12)
    ap.add_argument("--max-steps", type=int, default=None)
    args = ap.parse_args()
    r = run_b0(args.bin, args.out, args.tol, args.max_steps)
    print(f"[b0-lake] tier={r['tier']} max|eta-1|={r.get('max_eta_dev')} {r.get('reason','')}")
    return 0 if r["tier"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
