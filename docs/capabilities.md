# Frehg2 capabilities & the capability gate

This document is the authoritative list of what the realized Frehg2 solvers **actually
implement**. The production driver (`Orchestrator::initialize`) runs a **capability gate**
(`validateCapabilities` in `src/driver/Orchestrator.cpp`) that **fails fast** when a config
enables an option that is parsed by the schema but not implemented. This prevents the
"looks-wired-but-isn't" failure mode where a YAML key is read, stored, and then silently
ignored, so a run quietly computes physics you did not ask for.

If you hit `unsupported configuration: '<key>' ...`, the option is listed as *Not implemented*
below. Either disable/remove the key, or pick a supported alternative.

## Surface water (SWE)

Supported:
- Dynamic-wave semi-implicit SWE (legacy `difuwave == 0`).
- Manning friction, constant viscosity, waterfall correction, wet/dry handling.
- Rainfall (constant or time series) and **constant-rate** evaporation.
- Polygon boundary conditions over a region: `bc_discharge`, `bc_depth`, `bc_critical`.
- Polygon sources over a region: `inflow_rate`, `rainfall_rate`.

Not implemented (rejected by the gate when enabled):
- `surface_water.diffusive_wave: true` — only the dynamic wave is realized.
- `surface_water.wind.enabled: true` — wind forcing.
- `surface_water.subgrid.enabled: true` — subgrid bathymetry.
- `sources.surface.evaporation.model != 0` — evaporation **models** (only a constant/zero rate
  is applied).

## Groundwater (Richards / RE)

Supported:
- PCA predictor-corrector scheme (`groundwater.solver: pca`, `use_corrector: true`),
  vertical-only or full 3-D (`groundwater.full_3d: true`) with lateral Darcy flux.
- Per-cell heterogeneous soil (2-D per-column or 3-D per-cell `soil.map`).
- Boundary modes per face (`groundwater.bc_type_gw`, legacy indices `x-, x+, y-, y+, z+
  (bottom), z- (top)`):
  - **Lateral faces (x-, x+, y-, y+):** no-flux (`0`) only.
  - **Bottom (z+):** no-flux (`0`), Dirichlet (`1`, `groundwater.hbot`), fixed-flux (`2`,
    `groundwater.qbot`; `qbot > 0` = upward inflow into the domain).
  - **Top (z-):** no-flux (`0`), Dirichlet (`1`, `groundwater.htop`), fixed-flux (`2`,
    `groundwater.qtop` / partial-width `groundwater.recharge`; `qtop < 0` = downward
    infiltration). When coupled to surface water the top flux is the SW<->GW exchange.
- Partial-width top recharge over a polygon region: `groundwater.recharge: {rate, polygon}`.
- Extraction-well sink over a polygon region (`sources:` `extraction_well`), clamped at the
  per-cell residual moisture.

Not implemented (rejected by the gate when selected):
- **Lateral Dirichlet / fixed-flux** GW boundaries (`bc_type_gw[0..3] != 0`). The boundary
  value would be silently ignored, so these are rejected rather than mishandled.
- Iterative head solves: `groundwater.solver` other than `pca`, or `groundwater.iter_solve !=
  0` (Picard / Newton).

## Coupling, solute, density

Supported:
- Synchronous (sequential) and async (single-rank) SW<->GW coupling.
- Passive solute advection-diffusion-decay with rainfall mixing.
- **SW<->GW solute exchange**: dissolved solute is carried with the water that the coupling
  exchanges across the surface / top-GW-cell interface. The transfer is co-located with the
  water exchange (`Orchestrator::applyCouplingSoluteExchange` ->
  `applyInterfaceExchange`): infiltrating water carries the surface concentration into the top
  soil cell, seeping water carries the top-soil concentration into the surface water. Total
  solute mass is conserved to machine precision (gate: `test_solute_exchange`,
  `test_source_sink`).
- Solute **initial conditions** (`initial_conditions.solute.{surface,subsurface}`: constant /
  raster / formula / polygon) applied to the canonical conc storage at init.

Not implemented (rejected / restricted):
- `solute.baroclinic: true` — density-driven (baroclinic) coupling.
- Multi-rank async coupling falls back to synchronous (warned, numerically identical), pending
  the PETSc subcomm split.
- Solute is **passive** (no feedback on density / flow). Reactive multi-species chemistry is
  out of scope.

## Output

`output.variables` selects the written fields (validated against the enabled modules):
- Surface: `water_depth` (`depth`), `eta`, `u`/`uu` (velocity_u), `v`/`vv` (velocity_v).
- Subsurface: `head`, `moisture` (`water_content`), `qx`, `qy`, `qz`.
- `concentration` (`conc`) when solute is enabled.

An unknown name, or a field whose module is disabled, is a hard error. When `output.variables`
is absent, the default set is `water_depth, eta` (SW) and `head, moisture` (GW), plus
`concentration` when solute is on.
