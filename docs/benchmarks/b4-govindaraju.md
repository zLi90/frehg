# b4-govindaraju — overland-flow plane (Govindaraju)

**Suite:** SERGHEI test suite · **Registry gate:** `review` (tolerance: null) ·
**Frehg2 status:** ✅ runs (review-tier, Manning approximation)

## Physics

A rainfall-driven overland-flow plane: a 200 × 10 cell sloping DEM (cellsize 0.109725 m, the
plane drains in +x toward an outlet at x ≈ 21.9 m). A two-pulse hyetograph drives runoff; the
reference (`legacy/benchmarks/b4-govindaraju/ReferenceData/outflow.txt`) is the outlet discharge
hydrograph `Q(t)`. The legacy (SERGHEI) run uses **Chezy** friction (roughness 1.767), reflective
outer boundaries, an outlet polygon BC (`bctype 9`), and no infiltration (`InfModel: none`).

Hyetograph (converted to m/s in `benchmarks/b4-govindaraju/rain`): 50.8 mm/h baseline with
101.6 mm/h pulses over [601,1200] and [1801,2400] s, off after 2401 s.

## Frehg2 port (review-tier)

`benchmarks/b4-govindaraju/b4-govindaraju.yaml`:

- Surface-water only; DEM read as an ESRI raster (`dem.asc`).
- **Friction model approximation:** Frehg2 SWE uses **Manning** (n = 0.03 here), not Chezy. The
  hydrograph magnitude is therefore approximate; the registry gates this `review`.
- Outlet = polygon `bc_critical` (critical-depth weir) over the last downslope column.
- `min_depth = 1e-8` so the small per-step rain increment on the finely-celled plane accumulates
  rather than being clamped to dry.

`tests/integration/test_b4_govindaraju.cpp` runs a capped 600-step (30 s) window and asserts:
stable run; finite, non-negative depths; rainfall entered the system (`storage + outlet outflow
> 0`); and **no spurious water creation** (`storage + outflow ≤ rainfall delivered`). Observed:
storage ≈ 9.1e-3 m³, outflow ≈ 9.4e-5 m³, max depth ≈ 4.2e-4 m (= rain depth over 30 s).

## Three-tier acceptance

| Tier | Criterion |
|---|---|
| **PASS** | Outlet hydrograph matches `outflow.txt` to L2 < 1e-3 — needs a Chezy friction option (not implemented). |
| **REVIEW** *(current)* | Manning run is stable & mass-conserving; rain → runoff → outlet with the right qualitative hydrograph shape. |
| **FAIL** | Non-finite/negative depths, water created (`storage + outflow` exceeds rainfall), or no runoff. |
