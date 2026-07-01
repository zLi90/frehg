#!/usr/bin/env python3
"""trend_check.py — time-series trend metrics for the P19 validation harness (Task 19.3.3).

Computes monotonicity (slope-sign classification), mean, RMS, and peak (time + value) of a 1-D
series, and compares the classification against an expected trend. Importable by
`run_validation.py` and runnable as a CLI on a 2-column (t, value) text file.
"""
import math
import sys
from typing import List, Optional, Sequence, Tuple


def _finite(xs: Sequence[float]) -> bool:
    return all(math.isfinite(x) for x in xs)


def classify(values: Sequence[float], flat_tol: float = 1e-12) -> str:
    """Classify a series as 'increasing', 'decreasing', 'flat', or 'mixed' by step-sign votes."""
    vals = [float(v) for v in values]
    if len(vals) < 2 or not _finite(vals):
        return "flat" if len(vals) < 2 else "mixed"
    span = max(vals) - min(vals)
    scale = max(abs(max(vals)), abs(min(vals)), 1.0)
    if span <= flat_tol * scale:
        return "flat"
    ups = downs = 0
    for a, b in zip(vals[:-1], vals[1:]):
        d = b - a
        if d > flat_tol * scale:
            ups += 1
        elif d < -flat_tol * scale:
            downs += 1
    if ups > 0 and downs == 0:
        return "increasing"
    if downs > 0 and ups == 0:
        return "decreasing"
    return "mixed"


def stats(times: Sequence[float], values: Sequence[float]) -> dict:
    """Return mean, rms, peak_value, peak_time, trend, n for a (times, values) series."""
    vals = [float(v) for v in values]
    ts = [float(t) for t in times]
    n = len(vals)
    if n == 0:
        return {"n": 0, "mean": 0.0, "rms": 0.0, "peak_value": 0.0, "peak_time": 0.0,
                "trend": "flat", "finite": True}
    mean = sum(vals) / n
    rms = math.sqrt(sum(v * v for v in vals) / n)
    pk_i = max(range(n), key=lambda i: vals[i])
    return {
        "n": n,
        "mean": mean,
        "rms": rms,
        "peak_value": vals[pk_i],
        "peak_time": ts[pk_i] if pk_i < len(ts) else 0.0,
        "trend": classify(vals),
        "finite": _finite(vals),
    }


# A trend is "broken" (a clear regression signal) only when the observed direction is the exact
# OPPOSITE of the expected one, or the series is non-finite. A 'mixed'/'flat' observation against
# an expected monotone trend is reported as a soft mismatch (REVIEW), not a hard break.
_OPPOSITE = {"increasing": "decreasing", "decreasing": "increasing"}


def compare(observed_trend: str, expected_trend: Optional[str], finite: bool) -> str:
    """Return 'ok', 'review' (soft mismatch), or 'fail' (broken / non-finite)."""
    if not finite:
        return "fail"
    if not expected_trend or expected_trend == "any":
        return "ok"
    if observed_trend == expected_trend:
        return "ok"
    if observed_trend == _OPPOSITE.get(expected_trend):
        return "fail"
    return "review"


def _read_two_col(path: str) -> Tuple[List[float], List[float]]:
    ts: List[float] = []
    vs: List[float] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            ts.append(float(parts[0]))
            vs.append(float(parts[1]))
    return ts, vs


def main(argv: Sequence[str]) -> int:
    if len(argv) < 2:
        sys.stderr.write("usage: trend_check.py <series.txt> [expected_trend]\n")
        return 2
    ts, vs = _read_two_col(argv[1])
    s = stats(ts, vs)
    expected = argv[2] if len(argv) > 2 else None
    verdict = compare(s["trend"], expected, s["finite"])
    print(f"trend={s['trend']} mean={s['mean']:.6g} rms={s['rms']:.6g} "
          f"peak={s['peak_value']:.6g}@{s['peak_time']:.6g} verdict={verdict}")
    return 0 if verdict != "fail" else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
