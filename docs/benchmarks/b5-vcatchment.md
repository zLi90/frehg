# b5-vcatchment — V-catchment overland flow (Kollet et al.)

**Suite:** SERGHEI test suite · **Registry gate:** `review` (tolerance: null) ·
**Frehg2 status:** ✅ runs (review-tier, SW-only, uniform Manning)

## Physics

The classic V-catchment runoff test (Kollet et al., https://doi.org/10.1002/2016WR019191): a
101 × 55 real DEM (two hillslopes draining into a central channel that drains to an outlet), a
rainfall pulse (100 mm/h for the first 20 h of a 120 h simulation), and an outlet at the channel
mouth. The legacy references are an **inter-model intercomparison** (CATHY, HGS, ParFlow,
cast3m — `legacy/benchmarks/b5-vcatchment/discharge/*.csv`), not a single legacy snapshot, which
is why the registry gate is `review`.

## Frehg2 port (review-tier)

`benchmarks/b5-vcatchment/b5-vcatchment.yaml`:

- Surface-water only; DEM read as an ESRI raster (`dem.asc`).
- **Subsurface leg deferred:** the legacy case is a surface+subsurface coupling test
  (`async: 1`). The V-catchment is runoff-dominated, so SW-only captures the hydrograph shape for
  a review-tier check; full coupling on a 2-D catchment needs lateral Richards (see b3).
- **Uniform Manning** (n = 0.03); the legacy per-cell `roughness.input` field is not applied.
- Outlet = polygon `bc_critical` over the channel-outlet rows.

`tests/integration/test_b5_vcatchment.cpp` runs a capped 200-step (200 s) window and asserts:
correct grid (101 × 55); stable run; finite, non-negative depths; rainfall entered the system;
and no spurious water creation.

## Three-tier acceptance

| Tier | Criterion |
|---|---|
| **PASS** | Outlet hydrograph falls within the inter-model spread (CATHY/HGS/ParFlow) for the rain case — needs surface+subsurface coupling + per-cell roughness. |
| **REVIEW** *(current)* | SW-only run on the real DEM is stable & mass-conserving; rain → runoff → outlet. |
| **FAIL** | Non-finite/negative depths, water created, or the solver diverges on the real DEM. |
