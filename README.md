# Frehg2 User Manual

> Frehg2 is a C++17/Kokkos/MPI/PETSc upgrade of the Frehg surface-subsurface-solute
> model. It preserves Frehg's benchmark algorithms where required and adds selected
> SERGHEI-style capabilities such as asynchronous coupling, polygon-based
> boundary/source regions, ASCII raster inputs, non-uniform soil tables, flexible
> initial-condition utilities, and monitoring utilities.

## Contents

- [Model Capabilities](#model-capabilities)
- [Project Layout](#project-layout)
- [Build, Install, And Compile](#build-install-and-compile)
- [Run Frehg2](#run-frehg2)
- [Benchmark Validation](#benchmark-validation)
- [YAML Input Overview](#yaml-input-overview)
- [Complete YAML Field Reference](#complete-yaml-field-reference)
- [Initial Conditions](#initial-conditions)
- [Boundary Conditions](#boundary-conditions)
- [Polygon-Based BC And Source/Sink Regions](#polygon-based-bc-and-sourcesink-regions)
- [Spatial Input Files](#spatial-input-files)
- [Time-Series Input Files](#time-series-input-files)
- [Output Files](#output-files)
- [Known Current Limitations](#known-current-limitations)

## Model Capabilities

Frehg2 currently supports:

- Surface water: Frehg's semi-implicit shallow-water equation solver.
- Groundwater: Frehg's predictor-corrector Richards equation solver.
- Soil physics: uniform van Genuchten parameters plus optional non-uniform soil ID maps.
- Coupling utilities: Frehg-style synchronous helpers and SERGHEI-style asynchronous coupling.
- Boundary/source utilities: legacy boundary arrays plus optional polygon-defined regions.
- Initial-condition utilities: constant values, ESRI ASCII rasters, and HDF5 vector datasets.
- Monitoring utilities: point values, polygon summaries, mass-balance residuals, and HDF5 output.
- Output: HDF5 datasets and benchmark-compatible text output written beside `output.filename`.

The code is still release-blocked by benchmark comparison status. See
`docs/validation_report.md`.

## Project Layout

```text
frehg2/
├── CMakeLists.txt
├── README.md
├── upgrade_plan.md
├── benchmarks/
│   ├── b1-sw/
│   │   ├── b1-sw.yaml
│   │   ├── b1-input/
│   │   ├── out/
│   │   └── polygons/
│   └── b2-gw/
│       ├── b2-gw.yaml
│       ├── b2-input/
│       ├── out/
│       ├── polygons/
│       └── soil/
├── docs/
├── legacy/
│   ├── frehg.0/
│   ├── serghei/
│   └── benchmarks/
├── src/
│   ├── bc/
│   ├── core/
│   ├── coupling/
│   ├── io/
│   ├── re/
│   ├── solute/
│   └── swe/
└── tests/
```

## Build, Install, And Compile

### Prerequisites

This workspace expects the local dependency stack below. Do not replace these with
system packages unless the CMake configuration is intentionally updated.

- Compiler: `gcc-15` and `g++-15`
- CMake 3.16+
- MPICH: `/Users/zhili/Codes/local/mpich`
- Kokkos: `/Users/zhili/Codes/local/kokkos`
- PETSc: `/Users/zhili/Codes/local/petsc`
- HDF5: `/Users/zhili/Codes/local/hdf5`
- yaml-cpp: `/Users/zhili/Codes/local/yaml-cpp`

### Configure

```bash
cd /Users/zhili/Codes/frehg2
export CC=gcc-15
export CXX=g++-15
cmake -S . -B build
```

Common CMake options:

- `USE_KOKKOS=ON`: enable Kokkos parallel data structures.
- `USE_MPI=ON`: enable MPI support.
- `USE_PETSC=ON`: enable PETSc solver support.
- `USE_HDF5=ON`: enable HDF5 I/O support.
- `USE_YAML_CPP=ON`: enable yaml-cpp input parsing.
- `BUILD_TESTS=ON`: build unit and feature tests.
- `BUILD_DOCS=OFF`: build Doxygen documentation when Doxygen is available.

For example:

```bash
cmake -S . -B build -DBUILD_TESTS=ON -DBUILD_DOCS=OFF
```

### Compile

```bash
cmake --build build -j4
```

The executable is:

```text
build/src/frehg2
```

### Optional Install Step

The current project is primarily built and run from the build tree. If install
rules are added in the future, use:

```bash
cmake --install build --prefix /path/to/install
```

Until install rules are defined, run `build/src/frehg2` directly.

### Clean Rebuild

```bash
rm -rf build
export CC=gcc-15
export CXX=g++-15
cmake -S . -B build
cmake --build build -j4
```

## Run Frehg2

The executable currently accepts exactly one YAML configuration file:

```text
frehg2 <config.yaml>
```

The current `main.cpp` driver is a benchmark-oriented single-physics driver. It
runs the surface-water path when `surface_water.enable: true`; otherwise it runs
the groundwater path when `groundwater.enable: true`. Fully coupled production
runs require additional driver integration around the implemented coupling
utilities.

### Help

```bash
/Users/zhili/Codes/frehg2/build/src/frehg2 --help
```

### Surface-Water Benchmark

```bash
cd /Users/zhili/Codes/frehg2
./build/src/frehg2 benchmarks/b1-sw/b1-sw.yaml
```

### Groundwater Benchmark

```bash
cd /Users/zhili/Codes/frehg2
./build/src/frehg2 benchmarks/b2-gw/b2-gw.yaml
```

### MPI Execution

MPI support is part of the architecture. For an MPI-enabled build, use the MPI
launcher that matches the configured MPICH installation:

```bash
mpiexec -n 4 ./build/src/frehg2 benchmarks/b1-sw/b1-sw.yaml
```

Set `domain.use_mpi`, `domain.mpi_nx`, and `domain.mpi_ny` consistently with
the number of ranks. The current benchmark driver remains local/single-physics;
MPI execution support should be validated for new production driver paths.

## Benchmark Validation

Run all unit and feature tests:

```bash
ctest --test-dir build --output-on-failure
```

Run benchmark comparisons:

```bash
python3 tests/compare_benchmarks.py --benchmark b1-sw --new-output-dir benchmarks/b1-sw/out
python3 tests/compare_benchmarks.py --benchmark b2-gw --new-output-dir benchmarks/b2-gw/out
```

The documented validation target is relative L2 error below `1e-6` against
legacy benchmark references. Current release status and blockers are documented
in `docs/validation_report.md`.

## YAML Input Overview

Frehg2 reads YAML input files. Required top-level sections are:

```yaml
simulation: {}
domain: {}
time: {}
surface_water: {}
groundwater: {}
solute: {}
output: {}
monitor: {}
```

Optional Part 2 sections are:

```yaml
polygon_boundary_conditions: {}
polygon_sources: {}
```

Relative paths are resolved relative to the YAML file directory when using
`Config::baseDirectory()` helpers. The current benchmark driver also resolves
`simulation.legacy_input_dir` and `output.filename` relative to the YAML file.

### Minimal Surface-Water Example

```yaml
simulation:
  id: "surface-example"
  title: "Surface water example"
  author: "Frehg2"
  code_version: "2.0.0"
  legacy_input_dir: "b1-input/"

domain:
  nx: 1
  ny: 10
  nz: null
  dx: 80.0
  dy: 80.0
  dz: 0.1
  dz_incre: 1.0
  botZ: -3.0
  use_mpi: false
  mpi_nx: 1
  mpi_ny: 1
  bath_file: true
  actv_file: false

time:
  dt: 5.0
  Tend: 18000.0
  NT: 100
  dt_out: 1800
  Co_max: 2

surface_water:
  enable: true
  solver: "semi_implicit"
  init_eta: -2.0
  bc_type: [0, 0, 0, 0]
  rain_file: 1
  q_rain: 0.0
  q_evap: 0.0
  grav: 9.81
  min_depth: 1e-08
  manning: 0.019
  viscx: 1e-06
  viscy: 1e-06
  wtfh: 1e-08
  hD: 0.1

groundwater:
  enable: false
  bc_type: [0, 0, 0, 0, 0, 0]

solute:
  enable: false
  n_scalar: 0

output:
  format: "hdf5"
  filename: "out/output.h5"
  legacy_reference_dir: "out/"

monitor:
  n_monitor: 0
  points: []
```

### Minimal Groundwater Example

```yaml
simulation:
  id: "groundwater-example"
  legacy_input_dir: "b2-input/"

domain:
  nx: 1
  ny: 1
  nz: null
  dx: 1.0
  dy: 1.0
  dz: 0.01
  dz_incre: 1.0
  botZ: -1.0
  use_mpi: false
  mpi_nx: 1
  mpi_ny: 1
  bath_file: false
  actv_file: false

time:
  dt: 0.0001
  Tend: 46800
  NT: 100
  dt_out: 11700
  Co_max: 2

surface_water:
  enable: false
  bc_type: [0, 0, 0, 0]

groundwater:
  enable: true
  solver: "predictor_corrector"
  iter_solve: 0
  use_full3d: false
  dt_adjust: true
  follow_terrain: false
  sync_coupling: true
  async: true
  use_corrector: true
  post_allocate: false
  use_vg: true
  use_mvg: false
  aev: -0.02
  dt_min: 0.0001
  dt_max: 2.0
  Co_max: 2
  Ksx: 0.0
  Ksy: 0.0
  Ksz: 2.89e-06
  Ss: 1e-05
  soil_a: 1.43
  soil_n: 1.56
  wcs: 0.33
  wcr: 0.0
  init_wc: 0.033
  init_h: 1.0
  init_wt_rel: 0.5
  init_wt_abs: -0.75
  qtop: 0.0
  qbot: 0.0
  htop: 0.0
  hbot: 0.0
  qyp: 0.0
  qym: 0.0
  bc_type: [0, 0, 0, 0, 0, 1]

solute:
  enable: false
  n_scalar: 0

output:
  format: "hdf5"
  filename: "out/output.h5"

monitor:
  n_monitor: 0
  points: []
```

## Complete YAML Field Reference

This section lists known YAML fields used by the converted benchmark schema,
tests, utility APIs, and Part 2 additions. Fields marked as retained are kept
for compatibility with Frehg input mapping even if the current driver does not
fully consume them.

### `simulation`

- `simulation.id`: string model or benchmark identifier. Required.
- `simulation.title`: human-readable simulation title.
- `simulation.author`: author or creator string.
- `simulation.code_version`: Frehg2 version label.
- `simulation.legacy_input_dir`: directory containing converted legacy auxiliary files, such as
  `bath` and `rain`, relative to the YAML file.

### `domain`

- `domain.nx`: number of global x cells. Required.
- `domain.ny`: number of global y cells. Required.
- `domain.nz`: number of vertical cells. Optional; groundwater benchmark may compute it from
  `botZ`, `dz`, and `dz_incre` when set to `null`.
- `domain.dx`: x cell size. Required.
- `domain.dy`: y cell size. Required.
- `domain.dz`: nominal groundwater vertical spacing. Required.
- `domain.dz_incre`: vertical spacing multiplier. Required.
- `domain.botZ`: bottom elevation used by legacy groundwater grid construction.
- `domain.use_mpi`: enables MPI-domain behavior in compatible code paths.
- `domain.mpi_nx`: number of MPI partitions in x. Required by validation.
- `domain.mpi_ny`: number of MPI partitions in y. Required by validation.
- `domain.bath_file`: legacy bathymetry input flag. In current b1-sw driver, bathymetry is read
  from `simulation.legacy_input_dir/bath`.
- `domain.actv_file`: legacy active-cell input flag.

### `time`

- `time.dt`: main time step. Required.
- `time.Tend`: end time. Required.
- `time.NT`: retained legacy nominal time-step count.
- `time.dt_out`: output interval. Required.
- `time.Co_max`: Courant-like limiter retained from the legacy schema.

### `surface_water`

- `surface_water.enable`: enables surface-water driver path. Required.
- `surface_water.solver`: solver label; benchmark value is `"semi_implicit"`.
- `surface_water.difuwave`: retained diffusive-wave flag; diffusive wave is not part of the
  preserved Part 1 solver target.
- `surface_water.init_eta`: initial free-surface elevation.
- `surface_water.eta_file`: optional legacy free-surface initial-condition flag/path.
- `surface_water.uv_file`: optional legacy velocity initial-condition flag/path.
- `surface_water.bc_type`: four-entry legacy SW boundary code array in order
  `[x+, x-, y+, y-]`. Required.
- `surface_water.n_tide`: number of tide regions.
- `surface_water.tide_file`: tide data flag/path list.
- `surface_water.tide_dat_len`: tide record length list.
- `surface_water.tide_locX`: tide x index pairs.
- `surface_water.tide_locY`: tide y index pairs.
- `surface_water.init_tide`: initial tide values.
- `surface_water.rain_file`: rain file flag/path. Benchmark b1-sw uses legacy `1` and reads
  `simulation.legacy_input_dir/rain`.
- `surface_water.rain_dat_len`: rain record count.
- `surface_water.q_rain`: constant rain rate when file forcing is disabled or for polygon source
  compatibility.
- `surface_water.evap_file`: evaporation file flag/path.
- `surface_water.evap_model`: evaporation model selector retained from legacy input.
- `surface_water.q_evap`: constant evaporation rate.
- `surface_water.evap_dat_len`: evaporation record count.
- `surface_water.n_inflow`: number of inflow regions.
- `surface_water.inflow_file`: inflow data flag/path list.
- `surface_water.inflow_dat_len`: inflow record length list.
- `surface_water.inflow_locX`: inflow x index pairs.
- `surface_water.inflow_locY`: inflow y index pairs.
- `surface_water.init_inflow`: initial inflow values.
- `surface_water.grav`: gravitational acceleration.
- `surface_water.viscx`: x-direction viscosity.
- `surface_water.viscy`: y-direction viscosity.
- `surface_water.min_depth`: wet/dry minimum depth threshold.
- `surface_water.manning`: Manning roughness.
- `surface_water.wtfh`: wetting-front/head threshold used by the legacy solver.
- `surface_water.hD`: depth threshold parameter.
- `surface_water.rhoa`: air density.
- `surface_water.rhow`: water density.
- `surface_water.sim_wind`: wind forcing flag.
- `surface_water.wind_file`: wind forcing flag/path.
- `surface_water.wind_dat_len`: wind record count.
- `surface_water.init_windspd`: initial wind speed.
- `surface_water.init_winddir`: initial wind direction.
- `surface_water.Cw`: wind drag coefficient.
- `surface_water.CwT`: wind drag transition parameter.
- `surface_water.north_angle`: grid north angle.
- `surface_water.use_subgrid`: retained subgrid flag; subgrid is not implemented in Part 1.
- `surface_water.r_sub`: retained subgrid parameter.
- `surface_water.eta_sub_min`: retained subgrid parameter.
- `surface_water.eta_sub_max`: retained subgrid parameter.
- `surface_water.deta_sub`: retained subgrid increment.

### `groundwater`

- `groundwater.enable`: enables groundwater driver path. Required.
- `groundwater.solver`: solver label; benchmark value is `"predictor_corrector"`.
- `groundwater.iter_solve`: retained iterative solver selector. Part 1 uses `0`; Newton/Picard
  alternatives are not active targets.
- `groundwater.use_full3d`: enables 3D side-flow behavior.
- `groundwater.dt_adjust`: enables legacy adaptive groundwater stepping.
- `groundwater.follow_terrain`: terrain-following grid flag.
- `groundwater.sync_coupling`: retained synchronous-coupling selector.
- `groundwater.async`: Part 2 asynchronous-coupling option.
- `groundwater.use_corrector`: enables predictor-corrector correction.
- `groundwater.post_allocate`: retained legacy flag; disabled in Part 1 benchmark path.
- `groundwater.use_vg`: enables van Genuchten retention/conductivity relationships.
- `groundwater.use_mvg`: modified van Genuchten flag.
- `groundwater.aev`: air-entry value/model parameter.
- `groundwater.dt_min`: minimum groundwater time step.
- `groundwater.dt_max`: maximum groundwater time step.
- `groundwater.Co_max`: groundwater Courant-like adaptive limit.
- `groundwater.Ksx`: saturated hydraulic conductivity in x.
- `groundwater.Ksy`: saturated hydraulic conductivity in y.
- `groundwater.Ksz`: saturated hydraulic conductivity in z.
- `groundwater.Ss`: specific storage.
- `groundwater.soil_a`: uniform VG alpha parameter.
- `groundwater.soil_n`: uniform VG n parameter.
- `groundwater.wcs`: saturated water content.
- `groundwater.wcr`: residual water content.
- `groundwater.init_wc`: initial water content.
- `groundwater.init_h`: initial head mode/value.
- `groundwater.init_wt_rel`: relative water-table depth.
- `groundwater.init_wt_abs`: absolute water-table elevation.
- `groundwater.h_file`: optional legacy head initial-condition flag/path.
- `groundwater.wc_file`: optional legacy water-content initial-condition flag/path.
- `groundwater.qtop`: top flux value.
- `groundwater.qbot`: bottom flux value.
- `groundwater.htop`: top fixed-head value.
- `groundwater.hbot`: bottom fixed-head value.
- `groundwater.qyp`: y+ flux value.
- `groundwater.qym`: y- flux value.
- `groundwater.bc_type`: six-entry legacy GW boundary code array in order
  `[x+, x-, y+, y-, bottom(z+), top(z-)]`. Required.
- `groundwater.soil_map_file`: optional ESRI ASCII raster of surface soil IDs. The IDs are
  expanded through groundwater columns.
- `groundwater.soil_types`: optional list of non-uniform soil parameter records.
- `groundwater.soil_types[].id`: integer soil ID matched by `soil_map_file`.
- `groundwater.soil_types[].name`: optional soil name.
- `groundwater.soil_types[].soil_a`: VG alpha for this soil.
- `groundwater.soil_types[].soil_n`: VG n for this soil.
- `groundwater.soil_types[].wcs`: saturated water content for this soil.
- `groundwater.soil_types[].wcr`: residual water content for this soil.
- `groundwater.soil_types[].Ksx`: saturated conductivity x for this soil.
- `groundwater.soil_types[].Ksy`: saturated conductivity y for this soil.
- `groundwater.soil_types[].Ksz`: saturated conductivity z for this soil.

### `polygon_boundary_conditions`

- `polygon_boundary_conditions.enable`: enables polygon BC sections in utility/test paths.
- `polygon_boundary_conditions.surface`: list of surface polygon BC records.
- `polygon_boundary_conditions.groundwater`: list of groundwater polygon BC records.
- `polygon_boundary_conditions.<group>[].id`: user label.
- `polygon_boundary_conditions.<group>[].polygon_file`: polygon file path relative to YAML.
- `polygon_boundary_conditions.<group>[].type`: BC type string.
- `polygon_boundary_conditions.<group>[].value`: numeric BC value.

Recognized BC type strings in current tests/utilities:

- `fixed_water_level`
- `zero_gradient`
- `groundwater_head`

The C++ enum also defines:

- `fixed_flow_rate`
- `free_outflow`
- `tidal_water_level`
- `root_water_uptake`
- `gravity_drainage`

If additional strings are added to YAML, the loader that maps strings to
`BoundaryConditionType` must be updated accordingly.

### `polygon_sources`

- `polygon_sources.enable`: enables polygon source/sink sections in utility/test paths.
- `polygon_sources.surface`: list of surface source/sink records.
- `polygon_sources.groundwater`: reserved for groundwater source/sink records if used by a
  future driver path.
- `polygon_sources.<group>[].id`: user label.
- `polygon_sources.<group>[].polygon_file`: polygon file path relative to YAML.
- `polygon_sources.<group>[].type`: source/sink type string.
- `polygon_sources.<group>[].rate`: rate applied over the selected cells.

Recognized source/sink type strings are:

- `rainfall`: positive contribution.
- `inflow`: positive contribution.
- `evapotranspiration`: negative contribution.
- `pumping`: negative contribution.

### `solute`

- `solute.enable`: enables solute features where integrated.
- `solute.n_scalar`: number of scalar tracers.
- `solute.baroclinic`: retained baroclinic flag.
- `solute.superbee`: retained limiter flag.
- `solute.scalar_tide_file`: scalar tide file flag/path list.
- `solute.scalar_tide_datlen`: scalar tide data length list.
- `solute.scalar_inflow_file`: scalar inflow file flag/path list.
- `solute.scalar_inflow_datlen`: scalar inflow data length list.
- `solute.scalar_surf_file`: surface scalar IC/source file flag/path list.
- `solute.scalar_subs_file`: subsurface scalar IC/source file flag/path list.
- `solute.init_s_surf`: initial surface scalar values.
- `solute.init_s_subs`: initial subsurface scalar values.
- `solute.s_tide`: tide scalar values.
- `solute.s_inflow`: inflow scalar values.
- `solute.s_yp`: y+ scalar boundary values.
- `solute.s_ym`: y- scalar boundary values.
- `solute.difux`: scalar diffusivity x.
- `solute.difuy`: scalar diffusivity y.
- `solute.difuz`: scalar diffusivity z.
- `solute.disp_lon`: longitudinal dispersivity.
- `solute.disp_lat`: lateral dispersivity.

### `output`

- `output.format`: output format label. Current examples use `"hdf5"`.
- `output.filename`: HDF5 output file path. Required. Relative paths are resolved relative to the
  YAML file.
- `output.legacy_reference_dir`: legacy output directory metadata used for validation mapping.
- `output.variables`: list of requested variable names. Current benchmark driver writes its fixed
  benchmark variable set regardless of this list.

Common variable names:

- `water_depth`
- `water_surface_elevation`
- `velocity_x`
- `velocity_y`
- `hydraulic_head`
- `water_content`
- `darcy_flux_x`
- `darcy_flux_y`
- `darcy_flux_z`

### `monitor`

- `monitor.n_monitor`: number of monitor points. Required.
- `monitor.points`: list of point monitor records.
- `monitor.points[].name`: monitor label.
- `monitor.points[].type`: monitor type, such as `"cell"`.
- `monitor.points[].locX`: legacy x index.
- `monitor.points[].locY`: legacy y index.

## Initial Conditions

### Legacy Benchmark IC Fields

Surface-water ICs:

- `surface_water.init_eta`: constant initial free-surface elevation.
- `surface_water.eta_file`: optional legacy free-surface IC flag/path.
- `surface_water.uv_file`: optional legacy velocity IC flag/path.

Groundwater ICs:

- `groundwater.init_wc`: initial water content.
- `groundwater.init_h`: initial head mode/value.
- `groundwater.init_wt_rel`: relative water-table depth.
- `groundwater.init_wt_abs`: absolute water-table elevation.
- `groundwater.h_file`: optional legacy head IC flag/path.
- `groundwater.wc_file`: optional legacy water-content IC flag/path.

### Flexible IC Utility Sources

The `InitialCondition` utility supports three source types:

- `CONSTANT`: use one scalar value for every selected cell.
- `ASCII_RASTER`: read a 2D ESRI ASCII raster.
- `HDF5_VECTOR`: read a one-dimensional HDF5 dataset.

The internal utility spec is:

```cpp
InitialConditionSpec{
    variable,
    source,
    value,
    filename,
    dataset
}
```

For surface fields, raster size must match `domain.nx * domain.ny`. For
groundwater fields, a surface raster is expanded through all vertical cells in
each column.

## Boundary Conditions

### Legacy Surface-Water BC Array

```yaml
surface_water:
  bc_type: [0, 0, 0, 0]  # x+, x-, y+, y-
```

The array order follows legacy `bctype_SW`. Keep exactly four entries.

### Legacy Groundwater BC Array

```yaml
groundwater:
  bc_type: [0, 0, 0, 0, 0, 1]  # x+, x-, y+, y-, bottom(z+), top(z-)
```

Keep exactly six entries. For the b2-gw benchmark, top fixed head is represented
by the last entry and paired with `groundwater.htop`.

### Boundary Values

Groundwater scalar boundary values include:

- `groundwater.qtop`
- `groundwater.qbot`
- `groundwater.htop`
- `groundwater.hbot`
- `groundwater.qyp`
- `groundwater.qym`

Surface-water boundary and forcing values include:

- `surface_water.init_tide`
- `surface_water.init_inflow`
- `surface_water.q_rain`
- `surface_water.q_evap`

## Polygon-Based BC And Source/Sink Regions

Polygon sections are optional extensions that coexist with legacy fields.

### Surface Polygon BC Example

```yaml
polygon_boundary_conditions:
  enable: true
  surface:
    - id: "domain-zero-gradient"
      polygon_file: "polygons/domain.poly"
      type: "zero_gradient"
      value: 0.0
  groundwater: []
```

### Groundwater Polygon BC Example

```yaml
polygon_boundary_conditions:
  enable: true
  surface: []
  groundwater:
    - id: "top-head"
      polygon_file: "polygons/domain.poly"
      type: "groundwater_head"
      value: 0.0
```

### Polygon Source Example

```yaml
polygon_sources:
  enable: true
  surface:
    - id: "rain-zone"
      polygon_file: "polygons/domain.poly"
      type: "rainfall"
      rate: 1.0e-6
```

Positive source types are `rainfall` and `inflow`. Negative source types are
`evapotranspiration` and `pumping`.

### Polygon File Format

Polygon files contain one `x y` vertex per line. Lines may contain comments
after `#`. The polygon may be closed explicitly by repeating the first vertex.

```text
# Full domain polygon
0.0 0.0
80.0 0.0
80.0 800.0
0.0 800.0
0.0 0.0
```

The reader also supports an optional first row containing the vertex count:

```text
5
0.0 0.0
80.0 0.0
80.0 800.0
0.0 800.0
0.0 0.0
```

Cells are selected by cell center:

```text
x = (i + 0.5) * domain.dx
y = (j + 0.5) * domain.dy
```

For groundwater polygons, selected surface cells expand through all `k` cells in
the column.

## Spatial Input Files

### Legacy Auxiliary Files

The current benchmark driver reads legacy-style auxiliary files from
`simulation.legacy_input_dir`:

- `simulation.legacy_input_dir/bath`
- `simulation.legacy_input_dir/rain`

`bath` is a one-column surface-cell vector. It must contain `nx * ny` numeric
values:

```text
0.20
0.15
0.10
```

`rain` is read as a flat sequence of time/value pairs:

```text
0.0
5.5e-6
12000.0
5.5e-6
12001.0
0.0
18000.0
0.0
```

The rain forcing is linearly interpolated between pair records. This legacy
format is different from the generic `TimeSeries` reader described below.

### ESRI ASCII Raster Files

ASCII raster files are used by the raster reader, soil maps, and flexible IC
utilities.

```text
ncols 2
nrows 2
xllcorner 0.0
yllcorner 0.0
cellsize 1.0
NODATA_value -9999
3 4
1 2
```

Reader rules:

- Header keys must appear in order: `ncols`, `nrows`, `xllcorner` or `xllcenter`,
  `yllcorner` or `yllcenter`, `cellsize`, `NODATA_value`.
- `ncols`, `nrows`, and `cellsize` must be positive.
- Data rows are read from top row to bottom row, then stored internally with
  bottom row first.
- Values equal to `NODATA_value` are marked inactive.
- Active-cell IC and soil-map reads reject NODATA values in active model cells.

### Soil Map Example

```yaml
groundwater:
  soil_map_file: "soil/soil_map_uniform.asc"
  soil_types:
    - id: 0
      name: "legacy_uniform"
      soil_a: 1.43
      soil_n: 1.56
      wcs: 0.33
      wcr: 0.0
      Ksx: 0.0
      Ksy: 0.0
      Ksz: 2.89e-06
```

`soil_map_file` is a 2D surface raster of integer soil IDs. Frehg2 expands each
surface ID through its groundwater column.

## Time-Series Input Files

The `TimeSeries` reader accepts whitespace- or comma-delimited rows with two
columns:

```text
# time_seconds value
0, 0.0
3600, 1.0e-6
7200, 2.0e-6
```

or:

```text
0 0.0
3600 1.0e-6
7200 2.0e-6
```

Reader rules:

- Text after `#` is ignored.
- Times must be sorted ascending.
- Values are linearly interpolated between records.
- Query times before the first record use the first value.
- Query times after the last record use the last value.
- Empty files are rejected.

Legacy benchmark rainfall currently uses `simulation.legacy_input_dir/rain`
through the benchmark driver rather than the generic `TimeSeries` class.

## Output Files

### HDF5

`output.filename` controls the HDF5 file path:

```yaml
output:
  format: "hdf5"
  filename: "out/output.h5"
```

Relative output paths are resolved relative to the YAML file. The current driver
creates the parent directory if needed.

Surface-water datasets include paths like:

```text
/surface/water_depth/<time>
/surface/water_surface_elevation/<time>
/surface/velocity_x/<time>
/surface/velocity_y/<time>
/surface/inundation/<time>
```

Groundwater datasets include paths like:

```text
/groundwater/hydraulic_head/<time>
/groundwater/water_content/<time>
/groundwater/darcy_flux_x/<time>
/groundwater/darcy_flux_y/<time>
/groundwater/darcy_flux_z/<time>
/groundwater/zcell/0
```

### Benchmark-Compatible Text Output

The current benchmark driver also writes text files beside `output.filename`.
For `filename: "out/output.h5"`, text outputs are written under the same `out/`
directory.

Surface-water text files:

- `depth_<time>`
- `surf_<time>`
- `uu_<time>`
- `vv_<time>`
- `inun_<time>`

Groundwater text files:

- `head_<time>`
- `moisture_<time>`
- `qx_<time>`
- `qy_<time>`
- `qz_<time>`
- `zcell_0`
- `timestep`

## Known Current Limitations

- `v2.0.0` is not release-ready until benchmark comparison passes documented
  acceptance criteria.
- Some YAML fields are retained for legacy mapping and documentation but are not
  fully consumed by the current benchmark driver.
- The current executable is a single-physics benchmark driver: it chooses
  surface water or groundwater based on enable flags and does not yet run a full
  coupled SW+GW production simulation.
- Polygon BC/source sections are implemented and tested as utilities; full solver
  driver integration may require additional wiring for new production cases.
- Flexible IC utilities support constant, raster, and HDF5 vector sources, but
  YAML-to-IC driver integration should be extended for user-defined production
  runs.
- The native runtime format is YAML. Legacy `key=value` files should be converted
  rather than read directly.

## References

1. Li, S., Hodges, B.R., Shen, C. (2023). A 3D fully coupled surface-subsurface
   flow model with a robust and efficient semi-implicit solver. Journal of Hydrology.
2. Li, S., et al. Frehg reference materials in `legacy/`.
3. Caviedes-Voullième, D., et al. (2023). SERGHEI: A highly parallel resolved
   hydro-ecological model. Geoscientific Model Development.
4. Zheng, Y., et al. (2026). Coupled surface-subsurface asynchronous integration
   reference materials in `legacy/`.
