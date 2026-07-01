# b6-kuan — coupled surface-water / groundwater column (Kuan et al.)

**Suite:** legacy Frehg benchmark · **Registry gate:** `review` (tolerance: null) ·
**Frehg2 status:** ✅ runs (review-tier, PCA path, non-baroclinic)

> **Naming:** this benchmark is **`b6-kuan`** (the 1D–2D hybrid Kuan et al. test case), matching
> the legacy directory `legacy/benchmarks/b6-kuan`. (Earlier plan drafts used a different, wrong
> name; the correct identity is Kuan et al.)

## Physics

A coupled surface-water / groundwater column (1 × 68 surface cells over a 10-layer subsurface,
`dy = 0.05 m`, `dz = 0.04 m`, bottom at −0.4 m). Tidal/seepage forcing drives a coupled SW–GW
exchange with a salinity scalar. The legacy run
(`legacy/benchmarks/b6-kuan/input`, FrehdC key=value format) uses the **Newton** GW solver
(`iter_solve = 1`), **asynchronous** coupling (`sync_coupling = 0`, `n_substep = 50`), and a
**baroclinic** scalar (`n_scalar = 1`, `baroclinic = 1`). References:
`out-ss-syncV4` (steady) and `out-td-syncV4` (transient).

## Frehg2 port (review-tier)

`benchmarks/b6-kuan/b6-kuan.yaml`:

- Coupled SW + GW on the real 1 × 68 × 10 geometry (`bath` from `ex3_kuan_input/bath`), soil =
  the legacy VG parameters (α = 5.9, n = 2.68, Ks = 5.26e-3 m/s, θs = 0.46, θr = 0.04).
- **GW scheme approximation:** the production **PCA predictor-corrector** path is used
  (`solver: pca`, `use_corrector: true`). The legacy Newton solver (`iter_solve = 1`) is not
  implemented in Frehg2.
- **Solute/baroclinic disabled:** Frehg2 has a passive solute solver but no baroclinic SW–GW
  coupling; the salinity leg is deferred.
- **Coupling:** synchronous (sequential) for the review run; the single-rank async pipeline
  (`coupling.mode: async`) is available and numerically identical.

`tests/integration/test_b6_kuan.cpp` runs a capped 100-step (5 s) coupled window and asserts:
the run is coupled and stable; surface depths are finite & non-negative; groundwater
`wc ∈ [θr, θs]` everywhere; and the SW↔GW exchange volume is finite (the coupling is
conservative by construction, P6).

## Three-tier acceptance

| Tier | Criterion |
|---|---|
| **PASS** | Steady/transient head & seepage reproduce `out-ss-syncV4`/`out-td-syncV4` to L2 < 1e-2 — needs the Newton GW solver + baroclinic coupling. |
| **REVIEW** *(current)* | Coupled PCA run is stable, `wc ∈ [θr, θs]`, depths finite, conservative SW↔GW exchange. |
| **FAIL** | Divergence, out-of-range `wc`, non-finite depth, or a non-conservative exchange. |
