#!/usr/bin/env python3
"""run_validation.py — unified b1-b6 three-tier validation harness (P19).

The authoritative "is Frehg2 working?" check. For each of b1-sw, b2-gw, b3-kirkland,
b4-govindaraju, b5-vcatchment, b6-kuan it runs the production `frehg2` driver, then computes:

  * run health      (binary exit code),
  * finiteness      (all HDF5 field datasets finite),
  * water mass balance (closable for SW-only runs: final == initial + rain_in + inflow - outflow),
  * a trend signal  (surface/subsurface water-volume series direction vs trend_reference.yaml),
  * a primary metric per `validation_thresholds.yaml`:
      - legacy_l2_depth : relative L2 of surface depth vs a committed legacy reference (b1-sw),
      - review_physics  : registry `review` benchmarks — best tier is REVIEW (b2-gw, b3-b6).

and classifies each into PASS / REVIEW / FAIL. b0-lake is reported SEPARATELY (run_b0_lake.py).
Emits `validation_report.md` and `validation_report.json`. Exit code is non-zero if anything FAILs.

Usage:
    tools/run_validation.py [--bin build/frehg2] [--out build/validation] [--quick] [--only ID]
"""
import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time as _time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(REPO, "tools")
SUITE = ["b1-sw", "b2-gw", "b3-kirkland", "b4-govindaraju", "b5-vcatchment", "b6-kuan"]

# Tier severity (higher = worse) for combining sub-checks.
_RANK = {"PASS": 0, "REVIEW": 1, "FAIL": 2}
_NAME = {v: k for k, v in _RANK.items()}


def _worse(*tiers):
    return _NAME[max(_RANK[t] for t in tiers)]


def _abs(p):
    return p if os.path.isabs(p) else os.path.join(REPO, p)


def load_yaml(path):
    import yaml
    with open(path) as f:
        return yaml.safe_load(f)


def parse_summary(path):
    """Parse simulation_summary.txt into {key: value} (floats where possible)."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            key = parts[0]
            val = parts[1] if len(parts) > 1 else ""
            val = val.strip().strip('"')
            try:
                out[key] = float(val)
            except ValueError:
                out[key] = val
    return out


def cfg_path(bid):
    return os.path.join(REPO, "benchmarks", bid, bid + ".yaml")


def run_one(bid, entry, quick, quick_cap, binary, out_dir):
    """Run a single benchmark; return (rc, summary_path, h5_path, seconds)."""
    cfg = cfg_path(bid)
    doc = load_yaml(cfg)
    run_mode = entry.get("run", "capped")
    if run_mode == "full":
        doc.setdefault("time", {})["max_steps"] = 0  # uncap: let t_end govern (reach ref times)
    elif quick:
        doc.setdefault("time", {})["max_steps"] = int(quick_cap)
    else:
        doc.setdefault("time", {})["max_steps"] = int(entry.get("max_steps", 100))
    h5_path = os.path.join(out_dir, bid + ".h5")
    doc.setdefault("output", {})["filename"] = h5_path

    import yaml
    tmp = tempfile.NamedTemporaryFile("w", suffix=".yaml", dir=os.path.dirname(cfg), delete=False)
    yaml.safe_dump(doc, tmp)
    tmp.close()
    t0 = _time.time()
    try:
        rc = subprocess.run([_abs(binary), tmp.name], cwd=os.path.dirname(cfg)).returncode
    except FileNotFoundError:
        return 127, None, None, 0.0
    finally:
        os.unlink(tmp.name)
    summary_path = os.path.join(out_dir, "simulation_summary.txt")
    return rc, summary_path, h5_path, _time.time() - t0


def h5_finite_and_series(h5_path, dx, dy, dz, series_kind):
    """Return (finite, times, volumes). series_kind: surface_volume | subsurface_volume."""
    import h5py
    finite = True
    times, vols = [], []
    if not os.path.exists(h5_path):
        return False, [], []
    with h5py.File(h5_path, "r") as f:
        for top in ("/surface", "/subsurface"):
            grp = f.get(top)
            if grp is None:
                continue
            for fld in grp.keys():
                for tk in grp[fld].keys():
                    arr = grp[fld][tk][:]
                    if not all(math.isfinite(float(v)) for v in arr):
                        finite = False
        if series_kind == "surface_volume":
            grp = f.get("/surface/water_depth")
            factor = dx * dy
        else:
            grp = f.get("/subsurface/moisture")
            factor = dx * dy * dz
        if grp is not None:
            keys = sorted(grp.keys(), key=lambda k: float(k))
            for tk in keys:
                arr = grp[tk][:]
                times.append(float(tk))
                vols.append(float(sum(float(v) for v in arr)) * factor)
    return finite, times, vols


def legacy_l2_depth(h5_path, ref_path, dx, dy):
    """Relative L2 of /surface/water_depth vs a 'time v0 v1 ...' reference; None if unavailable."""
    import h5py
    if not (os.path.exists(h5_path) and os.path.exists(ref_path)):
        return None
    rows = []
    with open(ref_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split()
            rows.append((float(p[0]), [float(x) for x in p[1:]]))
    sse = ssr = 0.0
    used = 0
    with h5py.File(h5_path, "r") as f:
        grp = f.get("/surface/water_depth")
        if grp is None:
            return None
        for t, ref in rows:
            tk = str(int(t)) if float(t).is_integer() else None
            ds = grp.get(tk) if tk is not None else None
            if ds is None:
                continue
            model = [float(v) for v in ds[:]]
            if len(model) != len(ref):
                continue
            for m, r in zip(model, ref):
                sse += (m - r) ** 2
                ssr += r * r
            used += 1
    if used == 0 or ssr == 0.0:
        return None
    return math.sqrt(sse / ssr)


def mass_balance(summary, mass_fail):
    """Closable only for SW-only runs (no instrumented GW flux BC). Returns the relative residual
    and a `gate` that FAILs only when a *closable* budget exceeds `mass_fail` (so the SWE wet/dry
    conservation floor, ~1e-6, does not block a strict-L2 PASS). `band` is the informational
    PASS/REVIEW/FAIL classification per the plan's 1e-10 / mass_fail bands."""
    sw = summary.get("surface_water_present") in (1.0, "true", True)
    gw = summary.get("groundwater_present") in (1.0, "true", True)
    init = float(summary.get("water_initial_volume", 0.0) or 0.0)
    final = float(summary.get("water_final_volume", 0.0) or 0.0)
    rain = float(summary.get("water_rain_volume_in", 0.0) or 0.0)
    inflow = float(summary.get("water_polygon_inflow_volume", 0.0) or 0.0)
    outflow = float(summary.get("water_polygon_outflow_volume", 0.0) or 0.0)
    if not sw or gw:
        return {"closable": False, "rel": None, "gate": "PASS",
                "note": "GW flux BC not instrumented; budget informational"}
    expected = init + rain + inflow - outflow
    residual = final - expected
    denom = max(abs(init) + abs(rain) + abs(inflow), 1e-300)
    rel = abs(residual) / denom
    band = "PASS" if rel < 1e-10 else ("REVIEW" if rel < mass_fail else "FAIL")
    gate = "FAIL" if rel >= mass_fail else "PASS"
    return {"closable": True, "rel": rel, "residual": residual, "band": band, "gate": gate}


def evaluate(bid, entry, trend_ref, rc, summary_path, h5_path, dx, dy, dz):
    summary = parse_summary(summary_path)
    metric = entry.get("metric", "review_physics")
    series_kind = trend_ref.get(bid, {}).get("series", "surface_volume")
    expected_trend = trend_ref.get(bid, {}).get("expected_trend", "any")

    res = {"id": bid, "metric": metric, "rc": rc,
           "steps": summary.get("time_steps_completed"),
           "sim_time_s": summary.get("simulated_time_seconds")}

    if rc != 0:
        res.update({"tier": "FAIL", "reason": f"run rc={rc}"})
        return res

    finite, times, vols = h5_finite_and_series(h5_path, dx, dy, dz, series_kind)
    res["finite"] = finite

    import trend_check
    tstats = trend_check.stats(times, vols)
    tverdict = trend_check.compare(tstats["trend"], expected_trend, finite and tstats["finite"])
    res["trend"] = {"series": series_kind, "observed": tstats["trend"],
                    "expected": expected_trend, "verdict": tverdict,
                    "mean": tstats["mean"], "peak_value": tstats["peak_value"],
                    "peak_time": tstats["peak_time"]}
    trend_tier = {"ok": "PASS", "review": "REVIEW", "fail": "FAIL"}[tverdict]

    mb = mass_balance(summary, float(entry.get("mass_fail", 1e-2)))
    res["mass_balance"] = mb

    if not finite:
        res.update({"tier": "FAIL", "reason": "non-finite field values"})
        return res

    if metric == "legacy_l2_depth":
        ref = os.path.join(TOOLS, entry["reference"])
        l2 = legacy_l2_depth(h5_path, ref, dx, dy)
        res["legacy_l2"] = l2
        if l2 is None:
            metric_tier = "FAIL"
            res["reason"] = "could not compute legacy L2 (missing datasets/reference)"
        elif l2 < float(entry["strict"]):
            metric_tier = "PASS"
        elif l2 < float(entry["loose"]):
            metric_tier = "REVIEW"
        else:
            metric_tier = "FAIL"
        res["tier"] = _worse(metric_tier, mb["gate"], trend_tier)
    else:  # review_physics: registry-gated review; best attainable tier is REVIEW
        res["tier"] = _worse("REVIEW", mb["gate"], trend_tier)
        if entry.get("status") == "blocked":
            res["status"] = "blocked (full fidelity); reduced surrogate"
    return res


def write_reports(out_dir, results, b0, overall):
    os.makedirs(out_dir, exist_ok=True)
    payload = {"overall": overall, "b0_lake": b0, "benchmarks": results,
               "generated": _time.strftime("%Y-%m-%dT%H:%M:%S")}
    with open(os.path.join(out_dir, "validation_report.json"), "w") as f:
        json.dump(payload, f, indent=2)

    lines = ["# Frehg2 Unified Validation Report (P19)", "",
             f"Generated: {payload['generated']}", "",
             f"**Overall: {overall}**", "",
             "| Benchmark | Tier | Metric | Detail |", "|---|---|---|---|"]
    for r in results:
        if r["metric"] == "legacy_l2_depth":
            l2 = r.get("legacy_l2")
            detail = f"L2={l2:.3e}" if isinstance(l2, float) else "L2=n/a"
        else:
            mb = r.get("mass_balance", {})
            rel = mb.get("rel")
            detail = (f"mass_rel={rel:.2e}" if isinstance(rel, float) else "mass=informational")
        tr = r.get("trend", {})
        detail += f", trend={tr.get('observed','?')}/{tr.get('verdict','?')}"
        if r.get("status"):
            detail += f", {r['status']}"
        if r.get("reason"):
            detail += f", {r['reason']}"
        lines.append(f"| {r['id']} | {r['tier']} | {r['metric']} | {detail} |")
    lines += ["", "## b0-lake (separate well-balanced check)", "",
              f"- Tier: **{b0['tier']}** — max|eta-1| = {b0.get('max_eta_dev')} "
              f"(tol {b0.get('tol')})", ""]
    lines += ["## Tier meaning", "",
              "- **PASS**: strict legacy L2 met (b1-sw) / well-balanced to 1e-12 (b0-lake).",
              "- **REVIEW**: registry `review` benchmark ran stably with sound physics "
              "(mass balance, trends); human sign-off expected. Strict legacy parity for these "
              "is not attainable with the current solver scope (see docs/benchmarks/).",
              "- **FAIL**: crash, non-finite state, broken trend, or mass residual > threshold.", ""]
    with open(os.path.join(out_dir, "validation_report.md"), "w") as f:
        f.write("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default="build/frehg2")
    ap.add_argument("--out", default="build/validation")
    ap.add_argument("--quick", action="store_true", help="cap every benchmark to quick_max_steps")
    ap.add_argument("--only", default=None, help="run a single benchmark id")
    args = ap.parse_args()

    sys.path.insert(0, TOOLS)  # for `import trend_check` / `import run_b0_lake`
    th = load_yaml(os.path.join(TOOLS, "validation_thresholds.yaml"))
    trend_ref = load_yaml(os.path.join(TOOLS, "trend_reference.yaml")).get("benchmarks", {})
    entries = th["benchmarks"]
    quick_cap = int(th.get("quick_max_steps", 60))
    out_dir = _abs(args.out)
    os.makedirs(out_dir, exist_ok=True)

    ids = [args.only] if args.only else SUITE
    results = []
    for bid in ids:
        if bid not in entries:
            print(f"[validation] skip unknown {bid}")
            continue
        entry = entries[bid]
        cfg = load_yaml(cfg_path(bid))
        dom = cfg.get("domain", {})
        dx = float(dom.get("dx", 1.0)); dy = float(dom.get("dy", 1.0)); dz = float(dom.get("dz", 1.0))
        bench_out = os.path.join(out_dir, bid)
        os.makedirs(bench_out, exist_ok=True)
        print(f"[validation] running {bid} (metric={entry.get('metric')}) ...")
        rc, summ, h5, secs = run_one(bid, entry, args.quick, quick_cap, args.bin, bench_out)
        r = evaluate(bid, entry, trend_ref, rc, summ, h5, dx, dy, dz)
        r["seconds"] = round(secs, 2)
        print(f"[validation]   {bid}: {r['tier']}")
        results.append(r)

    # b0-lake (separate). Quick mode caps its steps too.
    import run_b0_lake
    b0 = run_b0_lake.run_b0(args.bin, os.path.join(out_dir, "b0-lake"),
                            tol=float(th.get("b0_lake", {}).get("strict", 1e-12)),
                            max_steps=quick_cap if args.quick else None)
    print(f"[validation] b0-lake: {b0['tier']} (max|eta-1|={b0.get('max_eta_dev')})")

    tiers = [r["tier"] for r in results] + [b0["tier"]]
    overall = "FAIL" if "FAIL" in tiers else ("REVIEW" if "REVIEW" in tiers else "PASS")
    write_reports(out_dir, results, b0, overall)
    print(f"[validation] OVERALL={overall}; reports -> {out_dir}/validation_report.{{md,json}}")
    return 1 if overall == "FAIL" else 0


if __name__ == "__main__":
    sys.exit(main())
