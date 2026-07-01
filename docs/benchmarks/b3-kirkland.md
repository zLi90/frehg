# b3-kirkland — 2-D layered-soil infiltration (Kirkland et al. 1992)

**Suite:** SERGHEI test suite · **Registry gate:** `review` (tolerance: null) ·
**Frehg2 status:** ✅ **FULL 2-D port** (3-D Richards + per-cell soil + fixed-flux top, PCA scheme)

## Physics

A 2-D (x–z) infiltration into a layered soil (Kirkland et al., 1992,
https://doi.org/10.1029/92WR00802). The domain is 50 × 30 cells (x × z, ny = 1), dx = dz = 0.1 m,
3 m deep. A partial-width fixed-flux source (`-5.787e-6 m/s`, downward) is applied at the top over
the central 30 columns; water infiltrates and spreads through a blocky two-class soil (a permeable
lens over a less-permeable band over a permeable base). Initial pressure head is a uniform `-500 m`
(very dry). The reference is the digitized Kirkland `h=0` and `h=-400` hydraulic-head contours from
the SERGHEI `makeplot.py` (`ref0` / `ref400`); the legacy `reference/` directory ships empty.

Legacy inputs: `legacy/benchmarks/b3-kirkland/input/` (`dem.input` 50×1, `subsurface.input`
`ndepth 30`, `vg.input` 2 classes, `soilID.input` 30×50 class grid, `gwbc.input` top fixed-flux
`-5.787e-6 m/s` over `polygontop.input`, `head.input` initial heads `-500`). The legacy `gw_scheme : 2`
is **Picard**; Frehg2 reproduces the case with the sanctioned **PCA predictor-corrector** instead.

## How Frehg2 represents it (previously blocked, now resolved)

Three capabilities were required; all three now exist:

1. **2-D lateral Richards** — the P23 3-D predictor-corrector (`groundwater.full_3d: true`) adds the
   lateral x Darcy flux to the corrector (the implicit head solve was already 3-D).
2. **Fully heterogeneous per-cell soil** — the P23 3-D soil map (`soil.map.layers`, one class raster
   per vertical layer) assigns different VG soils per (x, z) cell. Generated from the legacy
   `soilID.input` by `benchmarks/b3-kirkland/gen_soil_layers.py`.
3. **Partial-width top fixed-flux BC** — the legacy `bctype_GW[5]==2` / `qtop` path, implemented in
   `ReSolver` and wired as `groundwater.recharge` (a `rate` + `polygon` ring → per-column top flux).
   The predictor injects `qtop` (legacy `groundwater.c:666-669`); the PCA corrector's top-face flux
   is set to `qtop` (legacy's own commented `dqep`) so the non-iterative scheme conserves the
   prescribed recharge (Picard would re-solve to consistency; PCA needs the explicit top flux).

The **scheme is PCA, not Picard** (no iterative Picard/Newton, by policy). The match to the digitized
Kirkland contours is therefore approximate, not bit-exact.

## Result (committed run, t = 86 400 s)

* **Mass conservation:** groundwater storage increases by 0.1495 m³ vs. the prescribed recharge
  0.1500 m³ (`|rate| · 0.3 m² · 86 400 s`) — **99.6 %** conserved.
* **Contour match vs. digitized Kirkland reference:** the simulated `h=0` front matches `ref0` to
  **RMS ≈ 0.25 m**, and the `h=-400` front matches `ref400` to **RMS ≈ 0.29 m**, in a 3 m-deep
  domain. The 2-D plume forms under the source, spreads laterally, and is bounded by the soil
  layering — the expected qualitative behaviour.
* **State:** head and water content finite everywhere; `wc ∈ [θr, θs]` per cell.

Gated by `tests/integration/test_b3_kirkland.cpp` (full run, ~42 s): stability + bounded state,
mass conservation (within 5 %), a saturated plume under the source, and contour RMS < 0.5 m for both
levels.

## Three-tier acceptance

| Tier | Criterion |
|---|---|
| **PASS** | Reproduces the digitized Kirkland contours to a tight bit-exact band — would require Picard iteration (excluded by policy) and harmonic K-faces; not claimed. |
| **REVIEW** *(current)* | Full 2-D layered run is stable, mass-conservative (recharge to <5 %), `wc ∈ [θr, θs]`, forms the expected 2-D plume, and matches the Kirkland `h=0`/`h=-400` contours to ~0.25–0.29 m RMS (< 0.5 m). |
| **FAIL** | Diverges, non-finite / out-of-range `wc`, loses the recharge mass, or no plume / contour RMS ≥ 0.5 m. |
