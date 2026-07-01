# Frehg2 User Manual

Frehg2 is a production-grade, coupled **surface-water / groundwater / solute-transport**
numerical model. It solves the 2-D shallow-water equations (SWE) for overland/open-channel
flow, the 3-D Richards equation (RE) for variably-saturated subsurface flow, and a passive
advection–diffusion–decay equation for a dissolved solute, with mass-conservative coupling
between the surface and subsurface.

This manual is written for **new users**. It covers:

1. [What Frehg2 does (and does not) do](#1-what-frehg2-does)
2. [System requirements and dependencies](#2-system-requirements-and-dependencies)
3. [Building and compiling](#3-building-and-compiling)
4. [Running a simulation](#4-running-a-simulation)
5. [The configuration file — full field reference](#5-configuration-file-reference)
6. [Input data file formats](#6-input-data-file-formats)
7. [Complete sample input files](#7-sample-input-files)
8. [Output files and post-processing](#8-output-and-post-processing)
9. [Validation tools](#9-validation-tools)
10. [Capability limits (what gets rejected)](#10-capability-limits)
11. [Troubleshooting](#11-troubleshooting)

Related references: the authoritative generated schema is
[`docs/yaml_schema_v2.md`](yaml_schema_v2.md); supported/unsupported features are in
[`docs/capabilities.md`](capabilities.md); validation tiers are in
[`docs/validation.md`](validation.md).

---

## 1. What Frehg2 does

| Module | Equation | Enabled by | Notes |
|--------|----------|------------|-------|
| Surface water | 2-D semi-implicit shallow-water (dynamic wave) | `modules.surface_water: true` | Manning friction, rainfall/evaporation, wet/dry |
| Groundwater | 3-D Richards (predictor-corrector, "PCA") | `modules.groundwater: true` | van Genuchten / modified-VG soils, adaptive dt, per-cell soil map |
| Solute | passive advection–diffusion–decay | `modules.solute: true` + `solute.enabled: true` | upwind/MUSCL advection, implicit diffusion, rainfall mixing |
| Coupling | top-face Darcy exchange (SW↔GW) | both SW and GW enabled | mass-conservative water **and** solute exchange |

A simulation can use any subset of these modules. Common configurations:

- **Surface-water only** (`mode: surface_water`) — flooding, overland flow, channel routing.
- **Groundwater only** (`mode: groundwater`) — infiltration columns, recharge, water-table response.
- **Coupled** (`mode: coupled`) — integrated surface–subsurface hydrology.
- Add `modules.solute: true` to any of the above for transport of a dissolved constituent.

The single production executable is `frehg2`. It is driven entirely by **one YAML
configuration file** (schema version `2.0`) plus optional data files (DEM/bathymetry,
rainfall time series, soil-class raster).

---

## 2. System requirements and dependencies

Frehg2 is C++20 and links the following libraries:

| Dependency | Purpose | Notes |
|------------|---------|-------|
| **CMake** ≥ 3.20 | build system | |
| **C++20 compiler** | `g++-15` / `gcc-15` on this machine | must match the Kokkos C++ standard (20) |
| **Kokkos** | on-node parallelism (OpenMP/Serial; CUDA optional) | built with `Kokkos_CXX_STANDARD=20` |
| **MPI** | distributed-memory parallelism | MPICH on this machine; must match the MPI PETSc was built with |
| **PETSc** ≥ 3.25 | linear solvers (KSP) | the only sanctioned linear-solver backend |
| **HDF5** (C API, parallel) | output / checkpoints | `H5_HAVE_PARALLEL` |
| **yaml-cpp** | configuration parsing | |
| Python 3 + `h5py`, `PyYAML` (optional) | post-processing and validation harness | |

On the reference development machine (macOS / Darwin) all dependencies are installed under a
single prefix `/Users/zhili/Codes/local/`:

- Compilers: Homebrew `gcc-15` / `g++-15`.
- MPI: local MPICH (`/Users/zhili/Codes/local/bin/{mpicxx,mpicc,mpiexec}`), pinned to the
  build PETSc was compiled against. **Always use this `mpiexec` for multi-rank runs** —
  Homebrew's `mpiexec` runs the MPICH-linked binary as size-1 singletons.
- GPU: **not available on macOS** (CUDA-only). GPU code is compiled-out by default; real GPU
  execution is validated separately on a Linux/NVIDIA host.

---

## 3. Building and compiling

### 3.1 Standard (CPU/OpenMP) build

From the repository root:

```bash
cmake -S . -B build \
  -DCMAKE_C_COMPILER=gcc-15 \
  -DCMAKE_CXX_COMPILER=g++-15 \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces the driver executable `build/frehg2` plus the test suite. Configuration uses
the dependency prefix automatically; if your dependencies live elsewhere, pass
`-DCMAKE_PREFIX_PATH=/path/to/local` (and `-DHDF5_ROOT=...`).

`CMAKE_BUILD_TYPE` may be `Release` (optimized; recommended for production runs) or `Debug`
(assertions, no optimization; recommended while developing configs).

### 3.2 Run the test suite

```bash
ctest --test-dir build --output-on-failure          # all tests
ctest --test-dir build -R test_swe_b1 --output-on-failure   # one test by name
```

Tests are labelled `cpu`, `mpi`, `integration`, `regression`, `validation`, and `gpu`
(GPU tests are disabled unless built with CUDA). A clean checkout passes all enabled tests.

### 3.3 Optional build switches

| Option | Default | Effect |
|--------|---------|--------|
| `-DCMAKE_BUILD_TYPE=Release\|Debug` | `Release` (recommended) | optimization / assertions |
| `-DFREHG2_ENABLE_CUDA=ON` | `OFF` | enable the GPU backend and `gpu`-labelled tests (**Linux/NVIDIA only — fatal on macOS**) |
| `-DFREHG2_USE_LEGACY=ON` | `OFF` | build the archived legacy `frehg` code (opt-in; not part of Frehg2) |

### 3.4 Quick smoke test

```bash
./build/frehg2 --version
./build/frehg2 benchmarks/b1-sw/b1-sw.yaml
```

The second command runs the surface-water benchmark and prints, on success:

```
frehg2: run 'benchmarks/b1-sw/b1-sw.yaml' complete (NNN steps; summary -> .../simulation_summary.txt)
```

---

## 4. Running a simulation

### 4.1 Command line

```
frehg2 [options] config.yaml
```

| Option | Argument | Description |
|--------|----------|-------------|
| `-h`, `--help` | — | print usage and exit |
| `-v`, `--version` | — | print version and exit |
| `--restart FILE` | HDF5 checkpoint path | resume from a checkpoint instead of starting fresh |
| `--restart-time T` | seconds | simulation time stored in the checkpoint (defaults to 0) |

Any unrecognized `-...` flags (e.g. `-ksp_type gmres`, `-ksp_monitor`) are passed straight to
PETSc, so you can tune or debug the linear solver from the command line.

The **first non-flag argument** is the config file. Output paths inside the config are
resolved **relative to the directory containing the YAML file**, so runs are reproducible
regardless of your current working directory.

### 4.2 Serial run

```bash
./build/frehg2 path/to/config.yaml
```

"Serial" is simply one MPI rank. Output, summary, and monitor files are written under the
directory of `output.filename`.

### 4.3 Parallel (MPI) run

Use the **local MPICH** `mpiexec`:

```bash
/Users/zhili/Codes/local/bin/mpiexec -n 4 ./build/frehg2 path/to/config.yaml
```

The model owns its domain decomposition (the horizontal `x`/`y` plane is split across ranks;
the vertical `z` direction is always on-rank). The same executable and backend produce
bit-comparable results on 1, 2, and 4 ranks (rank-count equivalence is part of the test
suite). You may set an explicit process grid in the config:

```yaml
domain:
  mpi: { enabled: true, nx: 2, ny: 2 }   # 4 ranks arranged 2x2
```

If `domain.mpi.enabled` is false the decomposition is chosen automatically.

> **Note:** restart and the asynchronous coupling pipeline are **single-rank only**. A
> multi-rank async run automatically falls back to synchronous coupling (numerically
> identical) with a warning.

### 4.4 Restart / checkpointing

Enable periodic checkpoints in the config:

```yaml
time:
  dt_checkpoint: 3600      # write a checkpoint every 3600 s (0/null = off)
  max_checkpoints: 3       # keep only the 3 most recent checkpoint sets
```

Each checkpoint is an atomic, crash-safe HDF5 file (`<out>/...ckpt.<time>.h5`). Resume with:

```bash
./build/frehg2 config.yaml --restart out/run.ckpt.3600.h5 --restart-time 3600
```

A continuous run and a restarted run agree to machine precision. Restart is serial-only.

---

## 5. Configuration file reference

The configuration is a single YAML file. It **must** declare `schema_version: '2.0'` and
the five required top-level sections: `simulation`, `domain`, `time`, `modules`, `output`.
The driver rejects any other schema version or a missing required section.

Conventions in the tables below:

- **Units** use SI (m, s, m/s, etc.).
- "Default" means the value used when the key is omitted. Keys with no default are required
  in the context noted.
- All keys are lower-case; nested keys are shown as `parent.child`.

### 5.1 `schema_version` (required)

```yaml
schema_version: '2.0'
```

Must be the string `"2.0"`. To upgrade an old/experimental file, run
`./build/tools/migrate_yaml_v1_to_v2 in.yaml out.yaml`.

### 5.2 `simulation` (required)

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `id` | string | `frehg2` | run identifier; used in output/monitor/summary file names |
| `title` | string | — | free-text description (recorded in output provenance) |
| `mode` | string | `coupled` | `surface_water` \| `groundwater` \| `coupled` \| `solute`; warns/throws if it disagrees with `modules.*` |

### 5.3 `domain` (required)

Defines the grid geometry. The horizontal grid is `nx × ny` cells of size `dx × dy`. The
vertical grid (groundwater) has `nz` layers; the top layer is `dz` thick and successive
layers grow by `dz_incre`.

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `nx` | int | cells | — | interior cells in x (> 0) |
| `ny` | int | cells | — | interior cells in y (> 0) |
| `nz` | int | layers | — | vertical layers (> 0; used by groundwater) |
| `dx` | float | m | — | cell size in x |
| `dy` | float | m | — | cell size in y |
| `dz` | float | m | — | thickness of the top vertical layer |
| `dz_incre` | float | — | `1.0` | vertical layer growth factor (1.0 = uniform layers) |
| `botz` | float | m | `0.0` | domain bottom elevation (vertical datum) |
| `follow_terrain` | bool | — | `false` | terrain-following vertical grid |
| `x0`, `y0` | float | m | `0.0` | origin of the physical coordinate system (used by polygon regions) |
| `mpi.enabled` | bool | — | `false` | use an explicit MPI process grid |
| `mpi.nx`, `mpi.ny` | int | ranks | `1` | process-grid dimensions |
| `bathymetry.from_file` | bool | — | `false` | load bed elevation from a file (see §6.1) |
| `bathymetry.file` | string | path | `bath` | bathymetry/DEM file path (relative to YAML) |
| `bathymetry.format` | string | — | `auto` | `auto` \| `raster`/`esri`/`ascii_raster` \| `list` |

If `bathymetry.from_file` is false, the bed is flat at `botz`.

### 5.4 `time` (required)

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `dt` | float | s | — | base time step |
| `t_end` | float | s | — | simulation end time |
| `max_steps` | int | steps | `0` | hard cap on steps (0 = derive from `t_end`/`dt`) |
| `output_interval` | int | s | `0` | field-output and monitor cadence (0 = no field output) |
| `dt_checkpoint` | int/null | s | `0`/null | checkpoint cadence (0/null = off) |
| `max_checkpoints` | int/null | count | null | number of checkpoint sets retained |

### 5.5 `modules` (required)

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `surface_water` | bool | `false` | enable the SWE solver |
| `groundwater` | bool | `false` | enable the Richards solver |
| `solute` | bool | `false` | enable solute transport (also requires `solute.enabled: true`) |

### 5.6 `surface_water`

Used when `modules.surface_water` is true.

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `gravity` | float | m/s² | `9.81` | gravitational acceleration |
| `manning` | float | s·m^(−1/3) | `0.0` | Manning roughness coefficient |
| `min_depth` | float | m | `1e-8` | wet/dry threshold (below this a cell is dry) |
| `viscosity.x`, `viscosity.y` | float | m²/s | `0.0` | horizontal eddy viscosity |
| `h_diffusion_ref` | float | m | `0.1` | reference depth for horizontal diffusion |
| `waterfall_depth` | float | m | `1e-8` | waterfall (steep-front) correction threshold |
| `solver.ksp_type` | string | — | `cg` | PETSc KSP type for the SW solve |
| `solver.pc_type` | string | — | `jacobi` | PETSc preconditioner |
| `solver.rtol` | float | — | `1e-12` | KSP relative tolerance |

> `min_depth` matters on fine grids with slow rainfall: water below `min_depth` is treated as
> dry and cleared each step, so set it below the per-step rain increment for slow-rain cases.

The following keys are **parsed but rejected if enabled** (see §10): `diffusive_wave`,
`wind.enabled`, `subgrid.enabled`.

### 5.7 `groundwater`

Used when `modules.groundwater` is true. Requires `soil.types` (§5.10).

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `solver` | string | — | `pca` | only `pca` (predictor-corrector) is supported |
| `full_3d` | bool | — | `false` | 3-D Darcy assembly (lateral + vertical) vs vertical-only |
| `adaptive_dt` | bool | — | `true` | non-CFL adaptive time stepping |
| `use_corrector` | bool | — | `true` | flux-based θ corrector (the PCA scheme) |
| `use_vg` | bool | — | `true` | van Genuchten retention |
| `use_mvg` | bool | — | `false` | modified van Genuchten |
| `air_entry_value` | float | m | `-0.02` | MVG air-entry pressure head |
| `specific_storage` | float | 1/m | `1e-5` | specific storage Sₛ |
| `dt_min` | float | s | `1e-4` | adaptive-dt lower clamp |
| `dt_max` | float | s | `2.0` | adaptive-dt upper clamp |
| `co_max` | float | — | `2.0` | adaptive-dt change criterion |
| `htop`, `hbot` | float | m | `0.0` | Dirichlet head at top/bottom face (when BC mode = 1) |
| `qtop`, `qbot` | float | m/s | `0.0` | fixed flux at top/bottom face (when BC mode = 2; negative = downward) |
| `bc_type_gw` | int[6] | — | all `0` | per-face BC mode, order `[x-, x+, y-, y+, z+ bottom, z- top]` |
| `recharge.rate` | float | m/s | — | uniform/region top recharge flux (negative = downward infiltration) |
| `recharge.polygon` | list of `[x,y]` | m | — | restrict recharge to columns inside this ring (omit ⇒ whole top) |
| `solver.ksp_type` / `pc_type` / `rtol` | — | — | `cg` / `sor` / `1e-10` | PETSc solver controls |

**`bc_type_gw` modes:** `0` = no-flux, `1` = Dirichlet (uses `htop`/`hbot`), `2` = fixed-flux
(uses `qtop`/`qbot`). **Lateral faces (indices 0–3) support only no-flux (0)**; Dirichlet/
fixed-flux are accepted only on the bottom (index 4) and top (index 5) faces. Anything else is
rejected at startup.

### 5.8 `coupling`

Used when both surface water and groundwater are enabled.

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `mode` | string | — | `sequential` | `sync`/`sequential` (Gauss–Seidel) \| `async` (single-rank pipeline) |
| `surface_dt` | float | s | — | surface-water coupling window |
| `groundwater_dt` | float | s | — | groundwater coupling sub-step (GW catches up to the SW time) |

Coupling exchanges water **and** dissolved solute across the surface / top-soil interface,
conservatively, every coupling window.

### 5.9 `initial_conditions`

Each field accepts either a **scalar** (constant value) or a **mapping** with a `type`:

- `type: constant` + `value`
- `type: raster` (alias `file`) + `file` (+ `format`, `dataset` for HDF5)
- `type: formula` + `formula` (an expression in `x y z t`; supports `sin cos exp sqrt abs`)
- `type: polygon` + `default` + `values: [{name, value}, ...]` referencing
  `initial_conditions.regions`
- top-level `initial_conditions.type: restart` + `file` + `time` (resume from a checkpoint)

| Key | Default | Description |
|-----|---------|-------------|
| `surface_water.eta` | `domain.botz` | initial water-surface elevation (m) |
| `surface_water.u`, `surface_water.v` | `0.0` | initial velocities (m/s) |
| `groundwater.wc` | soil θ_r | initial volumetric water content |
| `groundwater.h` | `0.0` | initial pressure head (used instead of `wc` when head IC is raster/formula/polygon) |
| `solute.surface` | `0.0` | initial surface concentration (kg/m³) |
| `solute.subsurface` | `0.0` | initial subsurface concentration (kg/m³) |
| `regions` | — | list of `{name, vertices: [[x,y],...], value}` polygons for polygon ICs |

Legacy flags `*_from_file` (e.g. `surface_water.eta_from_file` + `eta_file`) map to raster ICs.

### 5.10 `soil` (required when groundwater is enabled)

`soil.types` is an **ordered list** of soil classes; the list index is the class id (class `0`
drives the uniform path when no per-cell map is given).

```yaml
soil:
  map:
    from_file: false        # set true to use a per-cell class raster
    file: soil_class.asc     # 2-D class-index file (list or ESRI raster)
    format: auto             # auto | raster/esri | list
    # layers: [c0.asc, c1.asc, ...]   # optional: one class file per vertical layer (3-D map)
  types:
    - id: 0
      theta_s: 0.40          # saturated water content
      theta_r: 0.08          # residual water content
      vg: { alpha: 1.0, n: 2.0 }       # van Genuchten parameters (1/m, -)
      k_sat: { x: 0.0, y: 1.16e-7, z: 1.16e-7 }   # saturated conductivity (m/s)
      # optional per-class overrides (else inherit groundwater.*):
      # specific_storage: 1.0e-5
      # use_vg: true
      # use_mvg: false
      # air_entry_value: -0.02
```

Required per class: `theta_s`, `theta_r`, `vg.{alpha,n}`, `k_sat.{x,y,z}`. The class-index
file values must be integers in `[0, number_of_types)`.

### 5.11 `solute`

Used when `modules.solute` is true. **`solute.enabled` is authoritative** — it must be true
for any solute work to happen (it overrides `modules.solute`).

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `enabled` | bool | — | `false` | master switch (must be true to run transport) |
| `c_rain` | float | kg/m³ | `0.0` | concentration of rainfall (rain mixing) |
| `k_decay` | float | 1/s | `0.0` | first-order decay constant |
| `D` | float | m²/s | `1e-9` | isotropic diffusion/dispersion coefficient |
| `advection_scheme` | string | — | `upwind` | `upwind` (1st order) \| `muscl` (2nd order) |
| `diffusion_scheme` | string | — | `implicit` | `implicit` (PETSc solve) \| `none` |
| `cfl_max` | float | — | `1.0` | advective CFL ceiling; a substep above this is refined |
| `min_depth` | float | m | `1e-8` | surface wet threshold for sources/advection |

> `solute.baroclinic: true` (density-driven coupling) is parsed but **rejected** — solute is
> passive only.

### 5.12 Polygon boundary conditions and sources

These are **sequences** at the top level (distinct from the `sources:` *map* used for uniform
rainfall/evaporation). Each entry has a `name`, a `type`, a `vertices` ring (≥ 3 `[x,y]`
points in physical coordinates), and a `rate`.

```yaml
boundaries:                       # polygon boundary conditions (surface)
  - name: outlet
    type: bc_critical             # bc_discharge | bc_depth | bc_critical
    vertices: [[990,0],[1000,0],[1000,100],[990,100]]
    rate: 0.0

sources:                          # polygon sources/sinks (SEQUENCE form)
  - name: inflow_patch
    type: inflow_rate             # inflow_rate | rainfall_rate | extraction_well
    vertices: [[0,0],[50,0],[50,50],[0,50]]
    rate: 0.5                      # m^3/s (inflow) / m/s (rainfall) / extraction (well)
```

| Boundary `type` | Meaning |
|------------------|---------|
| `bc_discharge` | prescribed discharge Q (m³/s); >0 outflow (clamped to available), <0 inflow |
| `bc_depth` | prescribed water depth/elevation override |
| `bc_critical` | critical-flow weir outlet |

| Source `type` | Meaning |
|----------------|---------|
| `inflow_rate` | volumetric inflow (m³/s) spread over the region |
| `rainfall_rate` | rainfall (m/s) applied only inside the region |
| `extraction_well` | subsurface water-content sink at the deepest GW cell (clamped at θ_r) |

Polygon regions and sources apply **every time step**, after the solve, on both serial and
parallel runs.

### 5.13 `monitoring` / `monitors`

Two equivalent forms are accepted. The benchmark-style `monitoring.points` list:

```yaml
monitoring:
  points:
    - { i: 0, j: 4 }          # surface probe at grid cell (i,j)
    - { i: 0, j: 4, k: 2 }    # subsurface probe (k => subsurface)
```

Or the richer `monitors` block with point probes and line-flux integrals:

```yaml
monitors:
  probes:
    - { name: gauge1, xyz: [400.0, 320.0], fields: [eta, depth, u, v] }
    - { name: cellA,  i: 10, j: 5, fields: [eta] }
  probes_subsurface:
    - { name: well1, xyz: [400.0, 320.0, -1.0], fields: [head, moisture] }
  lines:
    - { name: section, p0: [0,0], p1: [0,100], field: discharge }
```

Probe fields: surface `eta`, `depth`, `u`, `v`; subsurface `head`, `moisture`. Monitor output
is written at `time.output_interval` to a CSV (see §8.3). Monitoring requires
`output_interval > 0`.

### 5.14 `sources` (uniform forcing — map form)

When `sources:` is a **map** (not a sequence), it provides domain-wide rainfall and
evaporation:

| Key | Type | Unit | Default | Description |
|-----|------|------|---------|-------------|
| `surface.rainfall.from_file` | bool | — | `false` | read a `(time, rate)` time series from file |
| `surface.rainfall.file` | string | path | `rain` | rainfall time-series file (see §6.2) |
| `surface.rainfall.rate` | float | m/s | `0.0` | constant rainfall rate (used when `from_file` is false) |
| `surface.evaporation.rate` | float | m/s | `0.0` | constant evaporation rate |
| `surface.evaporation.model` | int | — | `0` | only `0` (constant rate) is supported |

### 5.15 `output` (required)

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `format` | string | `hdf5` | output backend (HDF5 is the default/only) |
| `filename` | string | — | output HDF5 path (relative to the YAML directory) |
| `io_mode` | string | `serial_gather` | `serial_gather` \| `parallel_collective` \| `file_per_rank` (see §8) |
| `variables` | string[] | sensible defaults | which fields to write (see below) |

If `output.filename` is omitted, no field output is written but the run still produces
`simulation_summary.txt` in the working directory.

**Recognized `output.variables`** (case-insensitive; aliases accepted):

| Canonical | Aliases | Module | Unit |
|-----------|---------|--------|------|
| `water_depth` | `depth` | surface | m |
| `eta` | `wse`, `water_surface_elevation` | surface | m |
| `u` | `uu`, `velocity_x`, `velocity_u` | surface | m/s |
| `v` | `vv`, `velocity_y`, `velocity_v` | surface | m/s |
| `head` | `pressure_head` | subsurface | m |
| `moisture` | `water_content`, `wc` | subsurface | m³/m³ |
| `qx`, `qy`, `qz` | — | subsurface | m/s (Darcy flux) |
| `concentration` | `conc` | both | kg/m³ |

Requesting a surface field with surface water disabled (or a subsurface field with groundwater
disabled) is an error. Omitting `variables` writes `water_depth`+`eta` (SW), `head`+`moisture`
(GW), and `concentration` when solute is on.

---

## 6. Input data file formats

All data file paths in the config are resolved **relative to the YAML file's directory**.

### 6.1 Bathymetry / DEM (`domain.bathymetry`)

Two formats, selected by `domain.bathymetry.format`:

- **`list`** — a plain text file of `nx*ny` numeric values in row-major order
  (`gi + gj*nx`, i.e. x fastest, then y). For a single-column domain (`nx == 1`) a file of
  `ny` values is also accepted. Values are bed elevations in metres.
- **`raster`/`esri`/`ascii_raster`** — an **ESRI ASCII raster** DEM with a header:

  ```
  ncols 200
  nrows 10
  xllcorner 0.0
  yllcorner 0.0
  cellsize 0.11
  nodata_value -9999
  <row 0 (north)>  ...
  <row nrows-1 (south)>  ...
  ```

  The raster must match `domain.nx`/`domain.ny`. ESRI rows run north→south; Frehg2 maps them
  so model `gj` increases northward. `format: auto` (the default) treats a `.asc` extension as
  a raster, otherwise a list.

> NODATA cells are currently filled with the lowest valid elevation and kept active (a warning
> is printed); true inactive-cell masking is a future capability.

### 6.2 Rainfall time series (`sources.surface.rainfall`)

When `from_file: true`, the file is a flat list of **alternating `(time, rate)` pairs**:

```
0       0.0
3600    1.0e-5
7200    2.0e-5
10800   0.0
```

Times are seconds; rates are m/s. The model linearly interpolates between samples. The file
must contain an even number of values.

### 6.3 Soil-class raster (`soil.map`)

A 2-D class-index file (list of `nx*ny` integers, or an ESRI `.asc` raster) whose values are
soil-class ids matching the order of `soil.types`. For a fully 3-D soil map, give one file per
vertical layer via `soil.map.layers: [layer0, layer1, ...]` (length `nz`).

### 6.4 Field-based initial conditions

`type: raster` IC fields use the same list / ESRI-raster formats as above; `format: hdf5` with
a `dataset` path reads from an HDF5 file.

---

## 7. Sample input files

### 7.1 Minimal surface-water run

```yaml
schema_version: '2.0'
simulation: { id: demo_sw, mode: surface_water }
domain:
  nx: 50
  ny: 50
  nz: 1
  dx: 10.0
  dy: 10.0
  dz: 0.1
  botz: 0.0
  bathymetry: { from_file: true, file: dem.asc, format: raster }
time: { dt: 1.0, t_end: 3600.0, max_steps: 0, output_interval: 300 }
modules: { surface_water: true, groundwater: false, solute: false }
surface_water: { gravity: 9.81, manning: 0.03, min_depth: 1.0e-6 }
initial_conditions:
  surface_water: { eta: 0.0 }
sources:
  surface:
    rainfall: { from_file: true, file: rain }
    evaporation: { rate: 0.0 }
boundaries:
  - name: outlet
    type: bc_critical
    vertices: [[490,0],[500,0],[500,500],[490,500]]
monitoring:
  points: [{ i: 25, j: 25 }]
output:
  format: hdf5
  filename: out/demo_sw.h5
  io_mode: serial_gather
  variables: [water_depth, eta, u, v]
```

### 7.2 Groundwater-only infiltration column

```yaml
schema_version: '2.0'
simulation: { id: demo_gw, mode: groundwater }
domain: { nx: 1, ny: 1, nz: 30, dx: 1.0, dy: 1.0, dz: 0.1, dz_incre: 1.0, botz: -3.0,
          follow_terrain: true }
time: { dt: 5.0, t_end: 46800.0, max_steps: 0, output_interval: 11700 }
modules: { surface_water: false, groundwater: true, solute: false }
groundwater:
  solver: pca
  full_3d: false
  adaptive_dt: true
  use_corrector: true
  use_vg: true
  specific_storage: 1.0e-5
  dt_min: 0.01
  dt_max: 20.0
  co_max: 2
  bc_type_gw: [0, 0, 0, 0, 0, 2]     # top face (index 5) = fixed flux
  qtop: -2.0e-7                       # downward infiltration (m/s)
soil:
  map: { from_file: false }
  types:
    - id: 0
      theta_s: 0.40
      theta_r: 0.08
      vg: { alpha: 1.0, n: 2.0 }
      k_sat: { x: 0.0, y: 1.16e-7, z: 1.16e-7 }
initial_conditions:
  groundwater: { wc: 0.10 }
output:
  format: hdf5
  filename: out/demo_gw.h5
  variables: [head, moisture, qz]
```

### 7.3 Coupled surface–subsurface with solute transport

```yaml
schema_version: '2.0'
simulation: { id: demo_coupled, mode: coupled }
domain: { nx: 20, ny: 20, nz: 10, dx: 5.0, dy: 5.0, dz: 0.05, dz_incre: 1.2, botz: -2.0 }
time: { dt: 2.0, t_end: 7200.0, max_steps: 0, output_interval: 600 }
modules: { surface_water: true, groundwater: true, solute: true }
surface_water: { gravity: 9.81, manning: 0.025, min_depth: 1.0e-6 }
groundwater:
  solver: pca
  full_3d: true
  adaptive_dt: true
  use_corrector: true
  use_vg: true
  specific_storage: 1.0e-5
  bc_type_gw: [0, 0, 0, 0, 0, 0]
coupling: { mode: sync, surface_dt: 2.0, groundwater_dt: 0.5 }
soil:
  map: { from_file: false }
  types:
    - id: 0
      theta_s: 0.42
      theta_r: 0.06
      vg: { alpha: 1.5, n: 1.8 }
      k_sat: { x: 1.0e-6, y: 1.0e-6, z: 1.0e-6 }
initial_conditions:
  surface_water: { eta: -0.5 }
  groundwater:   { wc: 0.20 }
  solute:        { surface: 10.0, subsurface: 0.0 }
sources:
  surface:
    rainfall: { from_file: false, rate: 1.0e-6 }
solute:
  enabled: true
  c_rain: 0.0
  k_decay: 0.0
  D: 1.0e-9
  advection_scheme: upwind
  diffusion_scheme: implicit
  cfl_max: 1.0
output:
  format: hdf5
  filename: out/demo_coupled.h5
  variables: [water_depth, eta, head, moisture, concentration]
```

Run any of these with:

```bash
./build/frehg2 demo_coupled.yaml
# or in parallel:
/Users/zhili/Codes/local/bin/mpiexec -n 4 ./build/frehg2 demo_coupled.yaml
```

The committed benchmark configs under `benchmarks/` (`b0-lake`, `b1-sw`, `b2-gw`,
`b3-kirkland`, `b4-govindaraju`, `b5-vcatchment`, `b6-kuan`) are full working examples.

---

## 8. Output and post-processing

### 8.1 HDF5 output structure

Field output goes to the single HDF5 file at `output.filename`. Datasets are organized as:

```
/simulation/metadata          # provenance attributes (git SHA, config hash, units, build info, ...)
/surface/<field>/<time>       # e.g. /surface/eta/600, /surface/water_depth/1200
/subsurface/<field>/<time>    # e.g. /subsurface/head/600, /subsurface/moisture/600
```

- `<field>` is the canonical name (`water_depth`, `eta`, `u`, `v`, `head`, `moisture`,
  `qx`/`qy`/`qz`, `concentration`).
- `<time>` is the snapshot time in seconds (whole seconds get an integer name like `600`;
  sub-second times get a trimmed-decimal name like `0.3`).
- Surface datasets are length `nx*ny` (global, row-major `gi + gj*nx`); subsurface datasets are
  `nx*ny*nz` (`gi + gj*nx + k*nx*ny`). On-disk data is always in **global physical order with
  no halo**, regardless of decomposition.
- Each field carries a `units` attribute. `/simulation/metadata` embeds `git_sha`,
  `config_sha256`, `build_type`, `compiler`, `kokkos_backend`, `solver_backend`, `mpi_ranks`,
  and `io_mode` for full reproducibility.

### 8.2 `io_mode` (write strategy)

| `output.io_mode` | Writer | Compression | Use when |
|------------------|--------|-------------|----------|
| `serial_gather` (default) | rank 0 gathers, one file | gzip-6 | small/medium runs; best storage |
| `parallel_collective` | all ranks, shared file, MPI-IO | none | large runs; best write bandwidth |
| `file_per_rank` | one file per rank | gzip-6 | very large scale; reassembled via stored index |

### 8.3 XDMF for ParaView / VisIt

On run completion, rank 0 writes an XDMF 3.0 sidecar (`<output>.xmf`) next to the HDF5 file.
Open the `.xmf` in **ParaView** or **VisIt** to view the surface fields as a time series —
no custom scripts needed.

### 8.4 Monitor CSV

If monitoring is configured and `output_interval > 0`, a CSV is written to
`<output_dir>/monitors/<simulation.id>.csv`. The header is
`time,<probe>.<field>,...,<line>.flux`; one row per output interval is flushed (and `fsync`-ed)
so the file survives a crash and can be appended on restart.

### 8.5 Simulation summary

Every run writes `<output_dir>/simulation_summary.txt` (or `./simulation_summary.txt` when no
output file is configured). It records step count, wall-clock timing per region
(`perf_*_seconds`), a water mass-balance budget (initial/final volume, rain in, polygon
inflow/outflow, coupling exchange), and — when solute is on — a solute mass balance
(`solute_initial_mass`, `solute_final_mass`, `solute_relative_mass_error`, etc.).

### 8.6 Reading results in Python

```python
import h5py
import numpy as np

nx, ny = 50, 50
with h5py.File("out/demo_sw.h5", "r") as f:
    # list snapshot times for a field
    times = sorted(f["/surface/eta"].keys(), key=float)
    print("eta snapshots:", times)

    # read one snapshot and reshape to the global grid (row-major gi + gj*nx)
    eta = f["/surface/eta/600"][:].reshape(ny, nx)

    # provenance / units
    md = f["/simulation/metadata"]
    print("git_sha:", md.attrs.get("git_sha"))
    print("eta units:", f["/surface/eta/600"].attrs.get("units"))

import matplotlib.pyplot as plt
plt.imshow(eta, origin="lower"); plt.colorbar(label="eta [m]"); plt.show()
```

Subsurface fields reshape to `(nz, ny, nx)`:

```python
with h5py.File("out/demo_gw.h5", "r") as f:
    head = f["/subsurface/head/11700"][:].reshape(30, 1, 1)  # (nz, ny, nx)
```

Monitor CSV:

```python
import pandas as pd
df = pd.read_csv("out/monitors/demo_sw.csv")
df.plot(x="time")
```

---

## 9. Validation tools

The unified harness runs all benchmarks through the production driver and classifies each as
PASS / REVIEW / FAIL:

```bash
python tools/run_validation.py --bin build/frehg2 --out build/validation   # all benchmarks
python tools/run_validation.py --quick                                     # fast smoke
python tools/run_b0_lake.py    --bin build/frehg2                          # lake-at-rest check
```

It writes `validation_report.{md,json}`. See [`docs/validation.md`](validation.md) for the
tiers and metrics. Other useful tools:

```bash
python tools/run_benchmark.py --list                       # show benchmarks + gates
python tools/perf_report.py --bin build/frehg2 --out build/perf   # performance sweep
./build/tools/migrate_yaml_v1_to_v2 old.yaml new.yaml      # upgrade an old config
python tools/gen_yaml_schema_doc.py                        # regenerate the schema reference
```

---

## 10. Capability limits

The driver runs a **capability gate** at startup and refuses configurations that enable a
feature it does not implement (rather than silently ignoring them). The error names the exact
key and points to [`docs/capabilities.md`](capabilities.md). Currently rejected when enabled:

- `surface_water.diffusive_wave` — only the dynamic-wave SWE is implemented.
- `surface_water.wind.enabled` — wind forcing is not implemented.
- `surface_water.subgrid.enabled` — subgrid bathymetry is not implemented.
- `sources.surface.evaporation.model` ≠ 0 — only constant-rate evaporation is implemented.
- `groundwater.solver` ≠ `pca` and `groundwater.iter_solve` ≠ 0 — Picard/Newton head solves
  are not implemented; use the PCA predictor-corrector scheme.
- `groundwater.bc_type_gw` lateral faces (indices 0–3) ≠ 0 — lateral GW Dirichlet/fixed-flux
  is not implemented; only no-flux on lateral faces.
- `solute.baroclinic: true` — density-driven coupling is not implemented (solute is passive).

Other current restrictions: restart and async coupling are single-rank; reactive/multi-species
chemistry is out of scope. GPU execution is compiled-out on macOS (CPU/OpenMP/MPI is the
validated path here).

---

## 11. Troubleshooting

| Symptom | Likely cause / fix |
|---------|--------------------|
| `unsupported configuration: '<key>' ...` | a capability-gate rejection — disable/remove the key (§10) |
| `schema_version` error or "missing required section" | add `schema_version: '2.0'` and all of `simulation/domain/time/modules/output`; or run `migrate_yaml_v1_to_v2` |
| `bathymetry file ... has N values; expected nx*ny` | the list file length must equal `nx*ny` (or `ny` when `nx==1`); for a DEM set `format: raster` |
| `raster DEM ... is AxB but domain is nx=..., ny=...` | the ESRI raster `ncols/nrows` must match `domain.nx/ny` |
| `soil class id ... out of range` | class-index file values must be `0..len(soil.types)-1` |
| `groundwater enabled but soil.types is empty` | add at least one `soil.types` entry |
| `KSPSolve` fails / singular system | check BCs leave the system well-posed (e.g. a closed box with `specific_storage: 0` and `nz=1` is singular — use a small positive `specific_storage`) |
| Multi-rank run behaves like N independent size-1 runs | you used the wrong `mpiexec`; use the local MPICH `mpiexec` |
| No HDF5 output produced | set `output.filename`, and `time.output_interval > 0` |
| Rainwater disappears on fine grids | lower `surface_water.min_depth` below the per-step rain increment |

For deeper architecture and developer topics, see
[`docs/architecture.md`](architecture.md) and [`docs/developer_guide.md`](developer_guide.md).
