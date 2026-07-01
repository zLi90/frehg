#!/usr/bin/env python3
"""Unit test for tools/trend_check.py (P19 Task 19.3.3). Pure Python (no h5py/frehg2 needed):
verifies the trend classifier and the expected-vs-observed comparison detect monotonicity breaks,
flats, and peaks. Asserts numerical bounds, not just "runs"."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "tools"))
import trend_check as tc  # noqa: E402


def main() -> int:
    # Classification.
    assert tc.classify([1, 2, 3, 4]) == "increasing"
    assert tc.classify([4, 3, 2, 1]) == "decreasing"
    assert tc.classify([1.0, 1.0, 1.0]) == "flat"
    assert tc.classify([1, 3, 2, 4]) == "mixed"

    # Stats: mean / rms / peak detection.
    s = tc.stats([0, 1, 2, 3], [0.0, 2.0, 4.0, 1.0])
    assert abs(s["mean"] - 1.75) < 1e-12, s["mean"]
    assert abs(s["peak_value"] - 4.0) < 1e-12
    assert abs(s["peak_time"] - 2.0) < 1e-12
    assert s["trend"] == "mixed"

    # compare(): expected vs observed. A clear OPPOSITE is a hard FAIL; a soft mismatch is REVIEW.
    assert tc.compare("increasing", "increasing", True) == "ok"
    assert tc.compare("decreasing", "increasing", True) == "fail"   # monotonicity break detected
    assert tc.compare("mixed", "increasing", True) == "review"      # soft mismatch
    assert tc.compare("increasing", "any", True) == "ok"
    assert tc.compare("increasing", "increasing", False) == "fail"  # non-finite -> fail

    # A NaN series is non-finite.
    s2 = tc.stats([0, 1], [float("nan"), 1.0])
    assert s2["finite"] is False

    print("trend_check unit test OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
