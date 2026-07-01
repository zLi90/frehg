# Frehg2 Unified Validation (P19)

`tools/run_validation.py` is the authoritative "is Frehg2 working?" check. It runs the production
`frehg2` driver on all six benchmarks and classifies each into one of three tiers, then emits
`validation_report.md` and `validation_report.json`. `b0-lake` (well-balanced lake-at-rest) is a
**separate** check, not part of the b1–b6 suite.

## Running

```bash
# Full suite (b1-sw runs to its real t_end; b2-b6 run to their per-benchmark capped step count)
python tools/run_validation.py --bin build/frehg2 --out build/validation

# Fast CI smoke (every benchmark capped to quick_max_steps; b1-sw still runs full for its L2 refs)
python tools/run_validation.py --bin build/frehg2 --out build/validation --quick

# A single benchmark
python tools/run_validation.py --only b4-govindaraju

# The separate well-balanced check (PASS criterion max|eta-1| <= 1e-12)
python tools/run_b0_lake.py --bin build/frehg2
```

The harness exits non-zero **only** on a `FAIL` tier, so it is safe to use directly as a CI gate.
It requires `h5py` and `PyYAML`.

## Tiers

| Tier | Meaning |
|---|---|
| **PASS** | Strict metric met: `b1-sw` relative depth L2 `< 1e-6` vs the committed legacy reference; `b0-lake` `eta` preserved to `1e-12`. |
| **REVIEW** | The benchmark ran stably with sound physics (finite fields, mass balance, expected trend direction), but cannot auto-PASS: the registry gates `b2-gw` and `b3-b6` as `review`, and the current solver scope cannot reproduce them to strict legacy L2 (see `docs/benchmarks/`). Maintainer sign-off expected. |
| **FAIL** | Crash / non-zero exit, non-finite field values, a broken trend (observed direction opposite the expected one), or a closable water-balance residual above the benchmark's `mass_fail`. |

Per-benchmark thresholds and run policy live in `tools/validation_thresholds.yaml`; trend
expectations in `tools/trend_reference.yaml`.

## What each metric checks

- **`legacy_l2_depth`** (`b1-sw`): relative L2 of `/surface/water_depth` against
  `benchmarks/b1-sw/frehg2_np1_depth_reference.txt` over the reference snapshot times. The physics
  matches the legacy-exact P4 port, so this is a true strict gate (`PASS < 1e-6`, `REVIEW < 1e-4`).
- **`review_physics`** (`b2-gw`, `b3-kirkland`, `b4-govindaraju`, `b5-vcatchment`, `b6-kuan`):
  run health + finiteness + trend + (where closable) water mass balance. Best attainable tier is
  `REVIEW`; a gross failure (crash, NaN/Inf, opposite trend, or mass blowup) drops it to `FAIL`.
- **`well_balanced`** (`b0-lake`, separate): `max|eta - 1| <= 1e-12` across all snapshots.

## Mass balance (Task 19.3.4)

The driver writes a water budget to `simulation_summary.txt`: `water_initial_volume`,
`water_final_volume`, `water_rain_volume_in`, `water_polygon_inflow_volume`,
`water_polygon_outflow_volume`, `water_polygon_well_volume`, `coupling_exchange_volume`, and the
`surface_water_present` / `groundwater_present` flags.

- **SW-only runs** have a *closable* budget: `final == initial + rain_in + inflow - outflow`. The
  harness reports the relative residual and bands it (`< 1e-10` PASS, `< mass_fail` REVIEW, else
  FAIL). The SWE semi-implicit wet/dry scheme has a small conservation floor (~1e-6 at rest,
  larger during initial dry-bed wetting — see `b5-vcatchment`'s tuned `mass_fail`), so SW mass
  balance normally lands in the REVIEW band, not at machine precision.
- **GW-present runs** (`b2-gw`, `b3`, `b6`) are *not* closable here: the groundwater top/bottom
  flux BC inflow is not instrumented in the summary, so the residual is reported as informational
  and does not gate the tier. (Element-wise GW legacy parity is enforced separately by ctest
  `test_re_b2_gw`.)

Failure-mode reading: large L2 + small mass error ⇒ likely IC; small L2 + large mass error ⇒
likely a flux bug; both large ⇒ structural.

## Trend checks (Task 19.3.3)

`tools/trend_check.py` builds a scalar time series (total surface water volume for SW runs, total
subsurface water volume for GW-only runs) from the HDF5 snapshots and classifies its direction
(`increasing` / `decreasing` / `flat` / `mixed`). A direction that is the exact opposite of the
expected one is a `FAIL`; a soft mismatch (`mixed`/`flat` vs a monotone expectation) is a `REVIEW`
note. It also reports mean, RMS, and peak time/value.

## CI integration (Task 19.3.5)

`.github/workflows/validation.yml` builds Frehg2, runs the unit/integration ctests, then the
unified harness. A `FAIL` fails the job (blocks merge); each `REVIEW` is surfaced as a GitHub
`::warning` the maintainer must acknowledge; the reports are uploaded as artifacts. The macOS/Linux
OpenMP runner covers all six benchmarks; a CUDA runner (for `gpu`-labeled benchmarks) is added once
available — none of b1–b6 are GPU-only today. Locally the same checks run as ctests
`test_validation_harness_quick`, `test_b0_lake_well_balanced`, and `test_trend_check`
(label `validation`).

## Current status (on the present code)

`b1-sw` **PASS** (L2 = 0 vs the committed reference), `b0-lake` **PASS** (`max|eta-1| ~ 2e-15`),
`b2-gw`/`b3-kirkland`/`b4-govindaraju`/`b5-vcatchment`/`b6-kuan` **REVIEW** (stable, sound physics).
Overall: **REVIEW** — expected, because five of six are registry-`review` benchmarks.
