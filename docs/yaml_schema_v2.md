# YAML Schema Reference (v2) — AUTHORITATIVE

> **Generated** by `tools/gen_yaml_schema_doc.py` from the committed example
> fixtures. Do not edit by hand; rerun the generator. This document is the full,
> authoritative reference (Appendix A of `INTEGRATED_PLAN.md` is a quick summary).

The production schema is `schema_version: "2.0"`. The single production driver
(`Orchestrator::initialize()`) refuses any config whose `schema_version` is not
`"2.0"` or that is missing a required top-level section. Field names are the P0-
frozen names; the v2 schema is additive only (no `grid.dims`, module-list, or
`io.dir`). Use `tools/migrate_yaml_v1_to_v2` to upgrade a v1/experimental YAML.

## Required top-level fields

- **`schema_version`** — production schema version; must be "2.0"
- **`simulation`**
- **`domain`**
- **`time`**
- **`modules`**
- **`output`**

Modules must include the three booleans `modules.{surface_water,groundwater,solute}`.

## Full key reference

Bold keys are required at the top level. `unit`/example/description are curated;
types and example values are read from the fixtures.

- **`schema_version`** (string, unit: —, e.g. `2.0`) — production schema version; must be "2.0"
- **`simulation`** (map)
  - `id` (string, unit: —, e.g. `b1-sw`) — run identifier (used for output/monitor file names)
  - `title` (string, unit: —, e.g. `Converted from legacy inp...`) — free-text description
  - `mode` (string, unit: —, e.g. `surface_water`) — surface_water | groundwater | coupled | solute
  - `legacy_finput` (string, unit: —, e.g. `b1-input/`)
- **`domain`** (map)
  - `nx` (int, unit: cells, e.g. `1`) — interior cells in x (>0)
  - `ny` (int, unit: cells, e.g. `10`) — interior cells in y (>0)
  - `nz` (int, unit: cells, e.g. `30`) — vertical layers (>0)
  - `dx` (float, unit: m, e.g. `80.0`) — cell size in x
  - `dy` (float, unit: m, e.g. `80.0`) — cell size in y
  - `dz` (float, unit: m, e.g. `0.1`) — top vertical layer thickness
  - `dz_incre` (float, unit: —, e.g. `1.0`) — vertical layer growth factor
  - `botz` (float, unit: m, e.g. `-3.0`) — domain bottom elevation (datum)
  - `follow_terrain` (bool, unit: bool, e.g. `true`) — terrain-following vertical grid
  - `mpi` (map)
    - `enabled` (bool, unit: bool, e.g. `false`) — use an explicit MPI process grid
    - `nx` (int, unit: ranks, e.g. `1`) — process-grid size in x
    - `ny` (int, unit: ranks, e.g. `1`) — process-grid size in y
  - `bathymetry` (map)
    - `from_file` (bool, unit: bool, e.g. `true`) — load bathymetry from a raster/list file
    - `source` (string, unit: —, e.g. `analytical`)
  - `active_mask` (map)
    - `from_file` (bool, unit: —, e.g. `false`)
- **`time`** (map)
  - `dt` (float, unit: s, e.g. `5.0`) — base time step
  - `t_end` (float, unit: s, e.g. `18000.0`) — simulation end time
  - `max_steps` (int, unit: steps, e.g. `100`) — step cap (0 = derive from t_end/dt)
  - `output_interval` (int, unit: s, e.g. `1800`) — field-output cadence (0 = off)
  - `dt_checkpoint` (null, unit: s, e.g. `null`) — checkpoint cadence (0/null = off)
  - `max_checkpoints` (null, unit: count, e.g. `null`) — retained checkpoint sets
- **`modules`** (map)
  - `surface_water` (bool, unit: bool, e.g. `true`) — enable the SWE solver
  - `groundwater` (bool, unit: bool, e.g. `false`) — enable the Richards solver
  - `solute` (bool, unit: bool, e.g. `false`) — enable solute transport
- `surface_water` (map)
  - `gravity` (float, unit: m/s^2, e.g. `9.81`) — gravitational acceleration
  - `manning` (float, unit: s/m^(1/3), e.g. `0.019`) — Manning roughness
  - `min_depth` (float, unit: m, e.g. `1e-08`) — wet/dry threshold
  - `viscosity` (map)
    - `x` (float, unit: —, e.g. `1e-06`)
    - `y` (float, unit: —, e.g. `1e-06`)
  - `waterfall_depth` (float, unit: m, e.g. `1e-08`) — waterfall correction threshold
  - `h_diffusion_ref` (float, unit: m, e.g. `0.1`) — reference depth for horizontal diffusion
  - `rho_air` (float, unit: —, e.g. `1.225`)
  - `rho_water` (float, unit: —, e.g. `998.0`)
  - `diffusive_wave` (bool, unit: —, e.g. `false`)
  - `bc_type_sw` (sequence, unit: —, e.g. `[4 item(s)]`)
  - `wind` (map)
    - `enabled` (bool, unit: —, e.g. `false`)
    - `from_file` (bool, unit: —, e.g. `false`)
    - `data_len` (int, unit: —, e.g. `100000`)
    - `init_speed` (float, unit: —, e.g. `0.0`)
    - `init_dir` (float, unit: —, e.g. `0.0`)
    - `cw` (float, unit: —, e.g. `0.0013`)
    - `cw_t` (float, unit: —, e.g. `5.0`)
    - `north_angle` (float, unit: —, e.g. `0.0`)
  - `subgrid` (map)
    - `enabled` (bool, unit: —, e.g. `false`)
    - `r_sub` (int, unit: —, e.g. `30`)
    - `eta_min` (float, unit: —, e.g. `0.0`)
    - `eta_max` (float, unit: —, e.g. `0.9`)
    - `deta` (float, unit: —, e.g. `0.05`)
  - `solver` (string, unit: —, e.g. `semi_implicit`)
- `groundwater` (map)
  - `solver` (string, unit: —, e.g. `pca`) — pca (predictor-corrector) is the production path
  - `full_3d` (bool, unit: bool, e.g. `true`) — 3-D Darcy assembly vs vertical-only
  - `adaptive_dt` (bool, unit: bool, e.g. `true`) — non-CFL adaptive time stepping
  - `use_corrector` (bool, unit: bool, e.g. `true`) — flux-based theta corrector (PCA)
  - `post_allocate` (bool, unit: —, e.g. `false`)
  - `use_vg` (bool, unit: bool, e.g. `true`) — van Genuchten retention
  - `use_mvg` (bool, unit: bool, e.g. `false`) — modified van Genuchten
  - `air_entry_value` (float, unit: m, e.g. `-0.02`) — MVG air-entry pressure head
  - `dt_max` (float, unit: s, e.g. `20.0`) — adaptive dt upper clamp
  - `dt_min` (float, unit: s, e.g. `0.01`) — adaptive dt lower clamp
  - `co_max` (int, unit: —, e.g. `2`) — adaptive-dt change criterion
  - `specific_storage` (float, unit: 1/m, e.g. `1e-05`) — specific storage Ss
  - `bc_type_gw` (sequence, unit: —, e.g. `[6 item(s)]`) — deprecated 6-int legacy BC [x+,x-,y+,y-,z+bottom,z-top]
- `coupling` (map)
  - `mode` (string, unit: —, e.g. `sync`) — sync | sequential | async
  - `surface_dt` (float, unit: s, e.g. `5.0`) — surface-water coupling window
  - `groundwater_dt` (float, unit: s, e.g. `0.01`) — groundwater coupling sub-step
- `initial_conditions` (map)
  - `surface_water` (map)
    - `eta` (float, unit: —, e.g. `-2.0`)
    - `eta_from_file` (bool, unit: —, e.g. `false`)
    - `uv_from_file` (bool, unit: —, e.g. `false`)
  - `groundwater` (map)
    - `wc` (float, unit: —, e.g. `0.5`)
    - `h` (float, unit: —, e.g. `1.0`)
    - `wt_rel` (float, unit: —, e.g. `1.0`)
    - `wt_abs` (float, unit: —, e.g. `-0.75`)
    - `h_from_file` (bool, unit: —, e.g. `false`)
    - `wc_from_file` (bool, unit: —, e.g. `false`)
  - `solute` (map)
    - `surface` (float, unit: —, e.g. `15.0`)
    - `subsurface` (float, unit: —, e.g. `35.0`)
- `boundary_conditions` (map)
  - `surface` (sequence, unit: —, e.g. `[0 item(s)]`)
  - `groundwater` (sequence, unit: —, e.g. `[6 item(s)]`)
  - `solute` (sequence, unit: —, e.g. `[0 item(s)]`)
- `sources` (map)
  - `surface` (map)
    - `rainfall` (map)
      - `from_file` (bool, unit: —, e.g. `true`)
      - `rate` (float, unit: —, e.g. `0.0`)
      - `data_len` (int, unit: —, e.g. `4`)
    - `evaporation` (map)
      - `from_file` (bool, unit: —, e.g. `false`)
      - `model` (int, unit: —, e.g. `0`)
      - `rate` (float, unit: —, e.g. `0.0`)
      - `data_len` (int, unit: —, e.g. `20000`)
    - `inflow` (map)
      - `count` (int, unit: —, e.g. `0`)
      - `loc_x` (sequence, unit: —, e.g. `[2 item(s)]`)
      - `loc_y` (sequence, unit: —, e.g. `[2 item(s)]`)
      - `from_file` (bool, unit: —, e.g. `false`)
      - `rate` (float, unit: —, e.g. `0.0`)
      - `data_len` (int, unit: —, e.g. `63`)
- `soil` (map)
  - `map` (map)
    - `from_file` (bool, unit: bool, e.g. `false`) — per-cell soil-class raster
  - `types` (sequence, unit: —, e.g. `[1 item(s)]`)
- `solute` (map)
  - `n_species` (int, unit: —, e.g. `0`)
  - `baroclinic` (bool, unit: —, e.g. `false`)
  - `advection_scheme` (string, unit: —, e.g. `upwind`)
  - `diffusion` (map)
    - `x` (float, unit: —, e.g. `1e-10`)
    - `y` (float, unit: —, e.g. `1e-10`)
    - `z` (float, unit: —, e.g. `1e-10`)
  - `dispersivity` (map)
    - `longitudinal` (float, unit: —, e.g. `0.1`)
    - `lateral` (float, unit: —, e.g. `0.01`)
  - `io` (map)
    - `tide_file` (bool, unit: —, e.g. `true`)
    - `tide_datlen` (int, unit: —, e.g. `732`)
    - `inflow_file` (bool, unit: —, e.g. `false`)
    - `inflow_datlen` (int, unit: —, e.g. `1`)
    - `surf_file` (bool, unit: —, e.g. `true`)
    - `subs_file` (bool, unit: —, e.g. `false`)
- `monitoring` (map)
  - `points` (sequence, unit: —, e.g. `[1 item(s)]`)
- **`output`** (map)
  - `format` (string, unit: —, e.g. `hdf5`) — hdf5 (default backend)
  - `filename` (string, unit: path, e.g. `out/b1-sw.h5`) — output file (relative to the YAML dir)
  - `variables` (sequence, unit: —, e.g. `[4 item(s)]`)
  - `legacy_foutput` (string, unit: —, e.g. `out/`)
- `validation` (map)
  - `reference_type` (string, unit: —, e.g. `legacy_text_snapshots`)
  - `tolerance` (float, unit: —, e.g. `1e-06`)
  - `variables` (sequence, unit: —, e.g. `[1 item(s)]`)
  - `well_balanced` (bool, unit: —, e.g. `true`)
  - `max_velocity` (float, unit: —, e.g. `1e-10`)
  - `max_eta_deviation` (float, unit: —, e.g. `1e-10`)
- `legacy_raw` (map)
  - `sim_id` (string, unit: —, e.g. `MAXP1`)
  - `finput` (string, unit: —, e.g. `b1-input/`)
  - `foutput` (string, unit: —, e.g. `out/`)
  - `NX` (int, unit: —, e.g. `1`)
  - `NY` (int, unit: —, e.g. `10`)
  - `botZ` (float, unit: —, e.g. `-3.0`)
  - `dx` (float, unit: —, e.g. `80.0`)
  - `dy` (float, unit: —, e.g. `80.0`)
  - `dz` (float, unit: —, e.g. `0.1`)
  - `dz_incre` (float, unit: —, e.g. `1.0`)
  - `use_mpi` (int, unit: —, e.g. `0`)
  - `mpi_nx` (int, unit: —, e.g. `1`)
  - `mpi_ny` (int, unit: —, e.g. `1`)
  - `dt` (float, unit: —, e.g. `5.0`)
  - `Tend` (float, unit: —, e.g. `18000.0`)
  - `NT` (int, unit: —, e.g. `100`)
  - `dt_out` (int, unit: —, e.g. `1800`)
  - `n_monitor` (int, unit: —, e.g. `1`)
  - `monitor_locX` (int, unit: —, e.g. `0`)
  - `monitor_locY` (int, unit: —, e.g. `4`)
  - `bath_file` (int, unit: —, e.g. `1`)
  - `grav` (float, unit: —, e.g. `9.81`)
  - `viscx` (float, unit: —, e.g. `1e-06`)
  - `viscy` (float, unit: —, e.g. `1e-06`)
  - `min_dept` (float, unit: —, e.g. `1e-08`)
  - `manning` (float, unit: —, e.g. `0.019`)
  - `wtfh` (float, unit: —, e.g. `1e-08`)
  - `hD` (float, unit: —, e.g. `0.1`)
  - `rhoa` (float, unit: —, e.g. `1.225`)
  - `rhow` (float, unit: —, e.g. `998.0`)
  - `sim_wind` (int, unit: —, e.g. `0`)
  - `wind_file` (int, unit: —, e.g. `0`)
  - `wind_dat_len` (int, unit: —, e.g. `100000`)
  - `init_windspd` (float, unit: —, e.g. `0.0`)
  - `init_winddir` (float, unit: —, e.g. `0.0`)
  - `Cw` (float, unit: —, e.g. `0.0013`)
  - `CwT` (float, unit: —, e.g. `5.0`)
  - `north_angle` (float, unit: —, e.g. `0.0`)
  - `sim_shallowwater` (int, unit: —, e.g. `1`)
  - `difuwave` (int, unit: —, e.g. `0`)
  - `init_eta` (float, unit: —, e.g. `-2.0`)
  - `eta_file` (int, unit: —, e.g. `0`)
  - `uv_file` (int, unit: —, e.g. `0`)
  - `bctype_SW` (sequence, unit: —, e.g. `[4 item(s)]`)
  - `n_tide` (int, unit: —, e.g. `0`)
  - `tide_file` (int, unit: —, e.g. `0`)
  - `tide_dat_len` (int, unit: —, e.g. `732`)
  - `tide_locX` (sequence, unit: —, e.g. `[2 item(s)]`)
  - `tide_locY` (sequence, unit: —, e.g. `[2 item(s)]`)
  - `init_tide` (float, unit: —, e.g. `-0.15`)
  - `evap_file` (int, unit: —, e.g. `0`)
  - `evap_model` (int, unit: —, e.g. `0`)
  - `q_evap` (float, unit: —, e.g. `0.0`)
  - `evap_dat_len` (int, unit: —, e.g. `20000`)
  - `rain_file` (int, unit: —, e.g. `1`)
  - `q_rain` (float, unit: —, e.g. `0.0`)
  - `rain_dat_len` (int, unit: —, e.g. `4`)
  - `inflow_locX` (sequence, unit: —, e.g. `[2 item(s)]`)
  - `inflow_locY` (sequence, unit: —, e.g. `[2 item(s)]`)
  - `n_inflow` (int, unit: —, e.g. `0`)
  - `inflow_file` (int, unit: —, e.g. `0`)
  - `init_inflow` (float, unit: —, e.g. `0.0`)
  - `inflow_dat_len` (int, unit: —, e.g. `63`)
  - `use_subgrid` (int, unit: —, e.g. `0`)
  - `r_sub` (int, unit: —, e.g. `30`)
  - `eta_sub_min` (float, unit: —, e.g. `0.0`)
  - `eta_sub_max` (float, unit: —, e.g. `0.9`)
  - `deta_sub` (float, unit: —, e.g. `0.05`)
  - `sim_groundwater` (int, unit: —, e.g. `0`)
  - `iter_solve` (int, unit: —, e.g. `0`)
  - `follow_terrain` (int, unit: —, e.g. `1`)
  - `sync_coupling` (int, unit: —, e.g. `1`)
  - `use_full3d` (int, unit: —, e.g. `1`)
  - `dt_adjust` (int, unit: —, e.g. `1`)
  - `use_corrector` (int, unit: —, e.g. `1`)
  - `post_allocate` (int, unit: —, e.g. `0`)
  - `use_vg` (int, unit: —, e.g. `1`)
  - `use_mvg` (int, unit: —, e.g. `0`)
  - `aev` (float, unit: —, e.g. `-0.02`)
  - `dt_max` (float, unit: —, e.g. `20.0`)
  - `dt_min` (float, unit: —, e.g. `0.01`)
  - `Co_max` (int, unit: —, e.g. `2`)
  - `Ksx` (float, unit: —, e.g. `0.0`)
  - `Ksy` (float, unit: —, e.g. `1.16e-07`)
  - `Ksz` (float, unit: —, e.g. `1.16e-07`)
  - `Ss` (float, unit: —, e.g. `1e-05`)
  - `soil_a` (float, unit: —, e.g. `1.0`)
  - `soil_n` (float, unit: —, e.g. `2.0`)
  - `wcs` (float, unit: —, e.g. `0.4`)
  - `wcr` (float, unit: —, e.g. `0.08`)
  - `init_wc` (float, unit: —, e.g. `0.5`)
  - `init_h` (float, unit: —, e.g. `1.0`)
  - `init_wt_rel` (float, unit: —, e.g. `1.0`)
  - `init_wt_abs` (float, unit: —, e.g. `-0.75`)
  - `h_file` (int, unit: —, e.g. `0`)
  - `wc_file` (int, unit: —, e.g. `0`)
  - `qtop` (float, unit: —, e.g. `0.0`)
  - `qbot` (float, unit: —, e.g. `0.0`)
  - `htop` (float, unit: —, e.g. `0.0`)
  - `hbot` (float, unit: —, e.g. `0.0`)
  - `qyp` (float, unit: —, e.g. `0.0`)
  - `qym` (float, unit: —, e.g. `0.0`)
  - `bctype_GW` (sequence, unit: —, e.g. `[6 item(s)]`)
  - `n_scalar` (int, unit: —, e.g. `0`)
  - `baroclinic` (int, unit: —, e.g. `0`)
  - `superbee` (int, unit: —, e.g. `0`)
  - `scalar_tide_file` (int, unit: —, e.g. `1`)
  - `scalar_tide_datlen` (int, unit: —, e.g. `732`)
  - `scalar_inflow_file` (int, unit: —, e.g. `0`)
  - `scalar_inflow_datlen` (int, unit: —, e.g. `1`)
  - `scalar_surf_file` (int, unit: —, e.g. `1`)
  - `scalar_subs_file` (int, unit: —, e.g. `0`)
  - `init_s_surf` (float, unit: —, e.g. `15.0`)
  - `init_s_subs` (float, unit: —, e.g. `35.0`)
  - `s_tide` (float, unit: —, e.g. `35.0`)
  - `s_inflow` (float, unit: —, e.g. `0.0`)
  - `s_yp` (float, unit: —, e.g. `0.0`)
  - `s_ym` (float, unit: —, e.g. `0.0`)
  - `difux` (float, unit: —, e.g. `1e-10`)
  - `difuy` (float, unit: —, e.g. `1e-10`)
  - `difuz` (float, unit: —, e.g. `1e-10`)
  - `disp_lon` (float, unit: —, e.g. `0.1`)
  - `disp_lat` (float, unit: —, e.g. `0.01`)
  - `actv_file` (int, unit: —, e.g. `0`)

## Source fixtures

This reference was generated from:
- `benchmarks/b1-sw/b1-sw.yaml`
- `benchmarks/b2-gw/b2-gw.yaml`
- `benchmarks/b0-lake/b0-lake.yaml`

## Migration (v1/experimental -> v2)

`tools/migrate_yaml_v1_to_v2 <in.yaml> [out.yaml]` applies the rename map (`src/io/YamlMigration.cpp`):

| v1 / experimental key | frozen v2 key |
|---|---|
| `time.max_step` | `time.max_steps` |
| `time.Tend` / top-level `Tend` | `time.t_end` |
| `time.dt_out` | `time.output_interval` |
| `grid.*` | `domain.*` |
| `domain.bot_z` | `domain.botz` |
| `output.directory` / `io.dir` | `output.filename` |
| BC/source `type: discharge\|depth\|critical` | `bc_discharge\|bc_depth\|bc_critical` |
| source `type: inflow\|rainfall\|well` | `inflow_rate\|rainfall_rate\|extraction_well` |
| `soil.uniform: true` (+ flat params) | `soil.types[0]` explicit class |

Archived v1 inputs live at `legacy/benchmarks/*/input.v1.yaml`.

