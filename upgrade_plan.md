# Frehg2 Upgrade Plan

> **Objective**: Upgrade the Frehg surface-subsurface-solute coupled numerical model from C/MPI/Makefile/LASPack to C++17/Kokkos/MPI/CMake/PETSc, adopt YAML input and HDF5 output, and add advanced features inspired by SERGHEI — while preserving Frehg's core Semi-Implicit SWE and Predictor-Corrector RE solver algorithms.

---

## Table of Contents

**Part 1: Architecture Upgrade — Same Algorithms, New Infrastructure**
1. [Phase 0: Legacy Behavior Audit](#1-phase-0-legacy-behavior-audit)
2. [Phase 1: Build System & Infrastructure](#2-phase-1-build-system--infrastructure)
3. [Phase 2: Core Data Structures](#3-phase-2-core-data-structures)
4. [Phase 3: I/O Layer — YAML Input & HDF5 Output](#4-phase-3-io-layer--yaml-input--hdf5-output)
5. [Phase 4: Surface Water Module — Semi-Implicit SWE](#5-phase-4-surface-water-module--semi-implicit-swe)
6. [Phase 5: Groundwater Module — Predictor-Corrector RE](#6-phase-5-groundwater-module--predictor-corrector-re)
7. [Phase 6: Surface-Subsurface Coupling — Frehg Algorithm](#7-phase-6-surface-subsurface-coupling--frehg-algorithm)
8. [Phase 7: Solute Transport Module — Frehg Algorithm](#8-phase-7-solute-transport-module--frehg-algorithm)

**Part 2: Algorithm Upgrade — New Features from SERGHEI**
9. [Phase 8: Asynchronous Coupling — from SERGHEI](#9-phase-8-asynchronous-coupling--from-serghei)
10. [Phase 9: Polygon-Based BC & Source/Sink — from SERGHEI](#10-phase-9-polygon-based-bc--sourcesink--from-serghei)
11. [Phase 10: Non-uniform Soil & Multi-VG — from SERGHEI](#11-phase-10-non-uniform-soil--multi-vg--from-serghei)
12. [Phase 11: Advanced Features — Flexible IC, Monitoring Upgrade](#12-phase-11-advanced-features--flexible-ic-monitoring-upgrade)
13. [Phase 12: Validation, Documentation & Release](#13-phase-12-validation-documentation--release)

14. [Appendix: File Format Specifications](#14-appendix-file-format-specifications)

---

## Design Principles

### Part 1 vs Part 2 — Clear Separation

| | Part 1: Architecture | Part 2: Algorithm |
|---|---|---|
| **Goal** | Port to new stack, preserve behavior | Add new capabilities |
| **Algorithms** | Identical to Frehg legacy | Changed/New (from SERGHEI) |
| **Validation** | Match legacy output byte-for-byte | Run benchmarks, verify correctness |
| **Input format** | YAML (1:1 mapping from legacy) | YAML (extended for new features) |

**Key rule**: Part 1 should produce results that are numerically identical (L2 error < 1e-6) to the legacy Frehg model for the same benchmark cases. If results differ, the *architecture port* has a bug.

### Algorithms to PRESERVE (from Frehg)

| Algorithm | Source | Notes |
|-----------|--------|-------|
| SWE Semi-Implicit Scheme | `frehg.0/src/shallowwater.c` | Keep core algorithm; remove Diffusive Wave and Subgrid Model |
| RE Predictor-Corrector | `frehg.0/src/groundwater.c` | Keep predictor-corrector; remove Picard, Newton, post-allocation |
| van Genuchten Model | `frehg.0/src/groundwater.c` | Keep VG and MVG models |
| Surface-Subsurface Coupling | `frehg.0/src/solve.c` | Keep Frehg's original coupling for Part 1 validation |
| Solute Transport | Frehg legacy | Keep Frehg's original solute algorithm for Part 1 validation |

### Algorithms to ADD (from SERGHEI — Part 2 only)

- Asynchronous coupling (rain → ponding → infiltration)
- Polygon-based flexible BC and source/sink application
- ASCII raster with NODATA support (ESRI format, for irregular domains)
- Non-uniform soil parameters (multi-VG, spatial soil type distribution via soil map)
- Flexible initial conditions (constant / file / multiple variables: h, wc, psi, saturation, etc.)

### New Infrastructure (both parts)

- C++17 + Kokkos for heterogeneous parallel computing (CPU/GPU portable)
- MPI for distributed memory parallelism
- PETSc for linear system solving (replace LASPack)
- YAML for input configuration (replace key=value format)
- HDF5 for Frehg2 output (legacy text outputs remain reference data for validation)
- CMake build system

### Development Environment

**Pre-installed dependencies** — Kokkos, MPICH, PETSc are installed at `/Users/zhili/Codes/local/`; **DO NOT** re-install or search system paths.

| Dependency | Path |
|---|---|
| Kokkos | `/Users/zhili/Codes/local/kokkos` |
| MPICH | `/Users/zhili/Codes/local/mpich` |
| PETSc | `/Users/zhili/Codes/local/petsc` |
| HDF5 | `/Users/zhili/Codes/local/hdf5` |
| yaml-cpp | `/Users/zhili/Codes/local/yaml-cpp` |

**Compiler**: `gcc-15` / `g++-15` (required):
```bash
export CC=gcc-15
export CXX=g++-15
cmake ..
```

---

## 1. Phase 0: Legacy Behavior Audit

> **Goal**: Fully document Frehg's current behavior before writing any new code. This audit is essential — it defines the "ground truth" that Part 1 must reproduce exactly.

### 0.1 Task: Document BC Code Meanings

**Sub-tasks:**

0.1.1. Read `frehg.0/src/shallowwater.c` and `frehg.0/src/groundwater.c` to catalog all BC type codes and their meanings.

0.1.2. Produce a BC code reference table documenting:
- Each integer BC code
- Its physical meaning (Dirichlet/Neumann/Freeflow/etc.)
- The corresponding variable being constrained
- How it is applied in the matrix system
- Where it appears in the code

0.1.3. Same for source/sink term codes.

**Deliverable**: `docs/legacy_audit/bc_code_reference.md`

---

### 0.2 Task: Document Index Conventions

**Sub-tasks:**

0.2.1. Trace how array indices work in Frehg:
- `i = 0..nx-1` (x-direction in cell)
- `j = 0..ny-1` (y-direction)
- Halo cell numbering
- Flat index: `idx = i + j * (nx + 2*halo)`
- 3D index for GW: `idx = i + j * nx + k * nx * ny` (verify)

0.2.2. Document MPI domain decomposition:
- How `mpi_nx` and `mpi_ny` divide the grid
- Local vs global indexing
- Halo exchange pattern

0.2.3. Document grid ordering:
- SW: 2D with halo width 1
- GW: 3D with halo width 1; vertical layers indexed bottom-to-top or top-to-bottom?

**Deliverable**: `docs/legacy_audit/index_conventions.md`

---

### 0.3 Task: Document State Variable Meanings

**Sub-tasks:**

0.3.1. Catalog every state variable array in Frehg:
- Name in code
- Physical meaning
- Units
- Initial value source
- Whether on host or kept with LASPack/Kokkos
- Time level (n, n-1, n+1)
- SW: `h`, `eta`, `hu`, `hv`, `z`, `roughness`
- GW: `h`, `wc`, `k`, `q`, `psi`

0.3.2. Document the flow of state variables through one time step for both SW and GW solvers.

**Deliverable**: `docs/legacy_audit/state_variables.md`

---

### 0.4 Task: Document Output Format (Legacy)

**Sub-tasks:**

0.4.1. For both b1-sw and b2-gw benchmarks, document:
- File naming convention (`depth_0`, `depth_1`, etc.)
- Legacy text output format (one floating-point value per line), dimensions, ordering
- What each file contains (which variable, at which time)
- How output time steps map to files

0.4.2. For b1-sw output (`legacy/benchmarks/b1-sw/out/`):
- `depth_*` files: 2D array of water depth, float/double, dimensions
- Any other output files

0.4.3. For b2-gw output (`legacy/benchmarks/b2-gw/out/`):
- `head_*`, `moisture_*`, `qz_*` files: 3D arrays, format details

0.4.4. Document monitoring output format.

**Deliverable**: `docs/legacy_audit/output_format.md`

---

### 0.5 Task: Map Every Legacy Input Field → YAML Schema

**Sub-tasks:**

0.5.1. Read benchmark inputs (`legacy/benchmarks/b1-sw/input`, `legacy/benchmarks/b2-gw/input`) line by line.

0.5.2. For every key=value pair, document:
- Key name
- Value type (int, real, string, array)
- Default value (if any)
- Rough meaning and where used in code
- Which module uses it (SW, GW, coupler)

0.5.3. Design YAML schema that captures every legacy input field 1:1:
- Section: `domain` (nx, ny, nz, dx, dy, dz, ...)
- Section: `time` (dt, Tend, dt_out, ...)
- Section: `surface_water` (all SW-related params)
- Section: `groundwater` (all GW-related params)
- Section: `solute` (all solute params if used)
- Section: `output` (output format settings)
- Section: `monitor` (monitoring point settings)

0.5.4. **Freeze the YAML schema after audit**. Any changes in Part 2 must extend (not break) this schema.

0.5.5. Preserve legacy positional BC ordering exactly. In particular, Frehg's groundwater
`bctype_GW` order is `[x+, x-, y+, y-, z+, z-]`, where `z+` is the bottom boundary
(`hbot`/`qbot`) and `z-` is the top/land-surface boundary (`htop`/`qtop`).

**Deliverable**: `docs/legacy_audit/yaml_schema.md` — complete 1:1 mapping table + frozen YAML schema

---

### 0.6 Task: Benchmark Conversion Rules

**Sub-tasks:**

0.6.1. Write a Python script `scripts/legacy_to_yaml.py` that converts old `input` files to new YAML format.

0.6.2. Run conversion for both benchmarks, producing `benchmarks/b1-sw/b1-sw.yaml` and `benchmarks/b2-gw/b2-gw.yaml`.

0.6.3. Verify manually that each field is correctly converted.

**Deliverable**:
- `scripts/legacy_to_yaml.py`
- Updated YAML files for both benchmarks

**Validation**: Visual review of converted YAML against legacy input.

---

## 2. Phase 1: Build System & Infrastructure

### 1.1 Task: Directory Structure

**Sub-tasks:**

1.1.1. Create directory tree:
```
frehg2/
├── benchmarks/
│   ├── b1-sw/           # SW rainfall-runoff
│   └── b2-gw/           # GW infiltration
├── cmake/                # CMake modules (FindKokkos, FindPETSc, etc.)
├── docs/
│   └── legacy_audit/    # Phase 0 audit documents
├── external/             # Third-party headers (if any)
├── scripts/              # Utility scripts (legacy_to_yaml.py, etc.)
├── src/
│   ├── core/             # Grid, Domain, State, types, MPI utils
│   ├── io/               # Config (YAML), Hdf5Writer, AsciiRaster, TimeSeries
│   ├── swe/              # Surface water module
│   ├── re/               # Richards equation module
│   ├── solute/           # Solute transport module
│   ├── coupling/         # Surface-subsurface coupling
│   └── bc/               # Boundary conditions (Part 2)
├── tests/
│   ├── core/
│   ├── io/
│   ├── swe/
│   ├── re/
│   ├── solute/
│   ├── coupling/
│   └── integration/
└── CMakeLists.txt
```

**Deliverables**: Complete directory structure, placeholder `.gitkeep` files.

**Unit Tests**: `tests/core/test_directory_structure.cpp` — verifies all directories exist.

---

### 1.2 Task: CMake Build System

**Sub-tasks:**

1.2.1. Write top-level `CMakeLists.txt`:
- C++17 standard
- Options: `USE_KOKKOS`, `USE_MPI`, `USE_PETSC`, `USE_HDF5`, `BUILD_TESTS`
- Find packages with explicit paths:
  ```cmake
  set(Kokkos_DIR /Users/zhili/Codes/local/kokkos)
  set(PETSC_DIR /Users/zhili/Codes/local/petsc)
  set(HDF5_ROOT /Users/zhili/Codes/local/hdf5)
  set(yaml-cpp_DIR /Users/zhili/Codes/local/yaml-cpp/share/cmake/yaml-cpp)
  ```
- Compiler: `gcc-15` / `g++-15`

1.2.2. Write `cmake/FindKokkos.cmake` and `cmake/FindPETSc.cmake` modules.

1.2.3. Write subdirectory `CMakeLists.txt` files (one per module with a stub library).

1.2.4. **Explicit discovery test**: After `find_package`, verify with:
  ```cmake
  # Verify Kokkos
  if(NOT Kokkos_FOUND)
    message(FATAL_ERROR "Kokkos not found at ${Kokkos_DIR}")
  endif()
  message(STATUS "Kokkos found: ${Kokkos_DIR}")
  # Same for PETSc, HDF5, MPI, yaml-cpp
  ```

1.2.5. Configure Kokkos with:
  ```cmake
  set(Kokkos_ENABLE_SERIAL ON)
  set(Kokkos_ENABLE_OPENMP ON)
  set(Kokkos_ENABLE_CUDA OFF)
  ```

**Deliverables**: Complete CMake system.

**Unit Tests**: `tests/core/test_cmake_config.cpp` — verifies all dependencies are discoverable and linkable.

**Validation**:
- `mkdir build && cd build && export CC=gcc-15 && export CXX=g++-15 && cmake ..` — must configure cleanly
- cmake output must explicitly confirm: "Kokkos found at /Users/zhili/Codes/local/kokkos", etc.
- `make -j4` — must compile empty target

---

### 1.3 Task: Core Type Definitions

**Sub-tasks:**

1.3.1. Create `src/core/define.hpp` — port from SERGHEI's `define.h`:
```cpp
#ifndef FREHG2_DEFINE_HPP
#define FREHG2_DEFINE_HPP

#include <Kokkos_Core.hpp>
#include <string>

namespace frehg2 {

// Real type
using real = double;

// Kokkos View types (host and device)
using realArr  = Kokkos::View<real*, Kokkos::HostSpace>;
using realArr2 = Kokkos::View<real**, Kokkos::HostSpace>;
using intArr   = Kokkos::View<int*, Kokkos::HostSpace>;
using realDeviceArr  = Kokkos::View<real*, Kokkos::DefaultExecutionSpace>;
using realDeviceArr2 = Kokkos::View<real**, Kokkos::DefaultExecutionSpace>;
using intDeviceArr   = Kokkos::View<int*, Kokkos::DefaultExecutionSpace>;

// Simulation mode flags
enum class SimMode { SW_ONLY, GW_ONLY, COUPLED, SOLUTE };

// BC types (extended during Part 2)
enum class BCType { DIRICHLET, NEUMANN, FREEFLOW, ZEROGRADIENT };

}  // namespace frehg2
#endif
```

1.3.2. Create `src/core/types.hpp` — domain and state structs.

**Deliverables**: `src/core/define.hpp`, `src/core/types.hpp`.

**Unit Tests**: `tests/core/test_types.cpp`.

**Validation**: Compile and run — must pass.

---

### 1.4 Task: Placeholder Main Program

**Sub-tasks:**

1.4.1. Create `src/main.cpp` — minimal main:
- Initialize Kokkos (`Kokkos::initialize()`)
- Initialize MPI (`MPI_Init()`)
- Initialize PETSc (`PetscInitialize()`)
- Parse command line: `./frehg2 config.yaml`
- `--help` flag

1.4.2. Create `src/core/Simulation.hpp/cpp` — stub `Simulation` class.

**Deliverables**: Compilable `main.cpp` with `--help`.

**Unit Tests**: `tests/core/test_main.cpp`.

**Validation**: `./frehg2 --help` — must print usage and exit cleanly.

---

## 3. Phase 2: Core Data Structures

### 2.1 Task: Grid and Domain Classes

**Sub-tasks:**

2.1.1. Create `src/core/Grid.hpp/cpp` — `Grid` class:
- Uniform grid: `nx, ny, nz, dx, dy, dz`
- Variable `dz` for GW (`dz_incre` multiplier)
- `nCell` (active cells), `nCellMem` (active + halo)
- Halo width = 1
- Methods: `getIndex(i,j,k)`, `isActive()`, `globalToLocal()`, `localToGlobal()`
- Index convention must match Phase 0 audit exactly

2.1.2. Create `src/core/Domain.hpp/cpp` — `Domain` class (2D surface):
- `z[i]` — bed elevation (Kokkos View)
- `area[i]` — cell area
- `actMask[i]` — active cell mask

2.1.3. Create `src/core/GwDomain.hpp/cpp` — `GwDomain` class (3D subsurface):
- `nz_glob`, `dz_multiplier`
- `soilID[i]` — soil type index (uniform in Part 1, extended in Part 2)

**Deliverables**: `Grid`, `Domain`, `GwDomain` classes.

**Unit Tests** (using gtest):
- `tests/core/test_grid.cpp`: Creates 10×10 grid, verifies `nx*ny == nCell`, tests halo cell indexing, tests `getIndex()` round-trip.
- `tests/core/test_domain.cpp`: Creates domain, verifies `actMask`, tests GW vertical layer creation.

**Validation**: All unit tests pass.

---

### 2.2 Task: State Variable Classes

**Sub-tasks:**

2.2.1. Create `src/core/State.hpp/cpp` — `State` class (surface water):
- `h[i]` — water depth
- `hu[i], hv[i]` — momentum
- `z[i]` — bed elevation
- `roughness[i]` — Manning's n
- `qss[i]` — surface-subsurface exchange flux

2.2.2. Create `src/core/GwState.hpp/cpp` — `GwState` class (groundwater):
- `h[i]` — hydraulic head (2 layers: old, new)
- `wc[i]` — water content (3 layers: old, new, excess)
- `k[i][3]` — hydraulic conductivity (x, y, z faces)
- `q[i][3]` — Darcy flux (x, y, z)

**Deliverables**: `State`, `GwState` classes with Kokkos Views.

**Unit Tests**:
- `tests/core/test_state.cpp`: Verifies View dimensions, tests deep_copy, tests state swap (old ↔ new).

**Validation**: All unit tests pass.

---

### 2.3 Task: MPI Communication Utilities

**Sub-tasks:**

2.3.1. Create `src/core/MpiComm.hpp/cpp`:
- Grid decomposition (`mpi_nx × mpi_ny`)
- Halo exchange for 2D and 3D arrays
- Gather/scatter for output
- Pattern must match Frehg's MPI behavior (from Phase 0 audit)

**Deliverables**: `MpiComm` class.

**Unit Tests**:
- `tests/core/test_mpi.cpp`: 4-process test on 10×10 grid, test halo exchange.

**Validation**: `mpirun -n 4 ./test_mpi` — must pass.

---

## 4. Phase 3: I/O Layer — YAML Input & HDF5 Output

### 3.1 Task: YAML Configuration Parser

**Sub-tasks:**

3.1.1. Create `src/io/Config.hpp/cpp` — `Config` class:
- Reads YAML using yaml-cpp
- Sections: `domain`, `time`, `surface_water`, `groundwater`, `solute`, `output`, `monitor`
- Methods: `get<T>(path)` with dot notation, e.g., `config->get<double>("domain.dx")`
- **Every field must map 1:1 to a legacy input field** (from Phase 0 audit)
- Schema validation: required fields check, type checking, default values

3.1.2. YAML schema is **frozen** based on Phase 0.5 output — all legacy fields are covered.

**Deliverables**: `Config` class.

**Unit Tests**:
- `tests/io/test_config.cpp`: Reads both benchmark YAMLs, verifies every legacy-mapped field is present and correctly typed.

**Validation**: Parse both benchmark YAML files without errors.

---

### 3.2 Task: HDF5 Output Writer

**Sub-tasks:**

3.2.1. Create `src/io/Hdf5Writer.hpp/cpp`:
- Creates HDF5 file with simulation metadata
- Writes time series datasets
- Parallel HDF5 with MPI
- Compression: gzip level 6
- Frehg2 writes simulation output in HDF5 format. Legacy text files are read only as
  validation references, not reproduced as the primary output format.

3.2.2. Dataset layout:
```
/simulation/
  /metadata (title, author, version, date)
  /domain (nx, ny, nz, dx, dy, dz, topography, active_mask)
  /time_series/timestamps
  /surface/water_depth, water_surface_elevation, velocity_x, velocity_y
  /subsurface/hydraulic_head, water_content, darcy_flux_x/y/z
  /monitor/point_001 (time series)
```

3.2.3. Create `src/io/Hdf5Reader.hpp/cpp` — for reading IC from HDF5.

**Deliverables**: `Hdf5Writer` and `Hdf5Reader`.

**Unit Tests**:
- `tests/io/test_hdf5.cpp`: Write then read, verify round-trip. Test parallel write (4 processes).

**Validation**: Write sample grid data, read with Python (h5py), compare.

---

### 3.3 Task: ASCII Raster Reader (NODATA Support)

**Sub-tasks:**

3.3.1. Create `src/io/AsciiRaster.hpp/cpp`:
- Read ESRI ASCII raster format (header: ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
- Store data in `realArr` (host)
- Set `actMask[i] = 0` where value == NODATA

3.3.2. Used for reading DEM, roughness, and (in Part 2) soil maps.

**Deliverables**: `AsciiRaster` reader.

**Unit Tests**:
- `tests/io/test_ascii_raster.cpp`: Read sample DEM with NODATA, verify header and values.

**Validation**: Read benchmark `bath` files, verify NODATA cells marked inactive.

---

### 3.4 Task: Time Series Reader

**Sub-tasks:**

3.4.1. Create `src/io/TimeSeries.hpp/cpp`:
- Read CSV/space-separated time series
- Linear interpolation: `getValueAt(time)`

**Deliverables**: `TimeSeries` reader.

**Unit Tests**:
- `tests/io/test_timeseries.cpp`: Read and interpolate.

**Validation**: Read benchmark rainfall data.

---

### 3.5 Task: Convert Legacy Inputs to YAML (One-Time)

**Sub-tasks:**

3.5.1. Run `scripts/legacy_to_yaml.py` (from Phase 0.6) to produce final YAML files.

3.5.2. Verify YAML files can be parsed by `Config` class.

3.5.3. **Note**: This is a one-time conversion for benchmark testing. Frehg2 does NOT need to support reading legacy `input` files natively.

**Deliverables**: Final `benchmarks/b1-sw/b1-sw.yaml` and `benchmarks/b2-gw/b2-gw.yaml`.

---

## 5. Phase 4: Surface Water Module — Semi-Implicit SWE

> **Critical**: This phase ports Frehg's Semi-Implicit SWE algorithm to C++/Kokkos/PETSc. The algorithm is **identical** to the legacy version — only the implementation language and linear solver are different.
>
> **Validation strategy**: Start with a **tiny 1D case** (10 cells) for rapid debugging, then progress to the full b1-sw benchmark.

### 4.0 Task: Tiny Deterministic Test — 1D SWE Matrix

**Sub-tasks:**

4.0.1. Create `tests/swe/test_1d_swe_matrix.cpp`:
- 10-cell 1D channel, flat bed, zero initial velocity, constant initial water depth
- Fixed-head BC at both ends
- Assemble the semi-implicit matrix (same algorithm as Frehg)
- Dump matrix to file
- Compare with matrix produced by: run legacy Frehg on same input, dump LASPack matrix
- **Every matrix element must match** (tolerance: 1e-12)

4.0.2. Create a helper script `scripts/dump_legacy_matrix.py` to extract LASPack matrix from legacy Frehg run.

**Deliverables**:
- `tests/swe/test_1d_swe_matrix.cpp`
- `scripts/dump_legacy_matrix.py`

**Validation**: Matrix element-by-element comparison passes → matrix assembly is correct.

---

### 4.1 Task: SWE Initialization

**Sub-tasks:**

4.1.1. Create `src/swe/SweSolver.hpp/cpp`:
- Preserve Frehg's Semi-Implicit Scheme algorithm from `shallowwater.c`
- **Do NOT implement** Diffusive Wave or Subgrid Model

4.1.2. Implement `initialize()`:
- Read IC from YAML: `init_eta` (constant or from file)
- Compute `h = max(0, eta - z)`, `hu = 0`, `hv = 0`
- BC handling: fixed `h`, fixed `hu/hv`, zero-gradient

**Deliverables**: `SweSolver` class with initialization.

**Unit Tests**:
- `tests/swe/test_swe_init.cpp`: Initialize with constant `eta`, verify `h = eta - z`.

**Validation**: Initialize b1-sw benchmark — compare initial `h` with reference.

---

### 4.2 Task: SWE Matrix Assembly

**Sub-tasks:**

4.2.1. Implement `assembleMatrix()` — Semi-Implicit discretization (port from `shallowwater_mat_coeff()`):
- Coefficient matrix `A` for `A * h_new = RHS`
- Use PETSc `Mat` for sparse matrix

4.2.2. Implement `assembleRHS()` — right-hand side (port from `build_shallowwater_system()`):
- Explicit terms from previous time step
- Source terms (rainfall, evaporation)

4.2.3. Index mapping: 2D grid → 1D PETSc matrix index must match Phase 0 audit.

**Deliverables**: Matrix assembly using PETSc.

**Unit Tests**:
- Re-run `tests/swe/test_1d_swe_matrix.cpp` — must still pass.
- `tests/swe/test_swe_matrix.cpp`: 5×5 grid, verify matrix values against manual calculation.

**Validation**: Dump matrix for b1-sw case, compare with legacy matrix dump.

---

### 4.3 Task: SWE Linear Solver (PETSc)

**Sub-tasks:**

4.3.1. Implement `solveMatrix()`:
- PETSc `KSP` with configurable solver/preconditioner
- Options: CG (+Jacobi, ILU), GMRES (+ILU)
- Solver tolerance: follow legacy LASPack accuracy (`SetRTCAccuracy(1e-8)`), e.g.
  `-ksp_rtol 1e-8` for PETSc runs unless Phase 0 audit shows a stricter equivalent is needed.

4.3.2. Implement `updateVelocity()` — explicit momentum update:
```
u_new = u_old - dt * (advection + pressure_gradient + friction + wind)
```

4.3.3. Implement legacy CFL diagnostics/checks:
```
dt_CFL = CFL * dx / max(|u| + sqrt(g*h))
```
- In Part 1, reproduce Frehg's existing surface-water time loop and do **not** introduce new
  SW adaptive time stepping before b1-sw matches legacy.
- Any SW adaptive time-step extension must be a separate Part 2/post-validation change.

**Deliverables**: PETSc-based linear solver, velocity update, legacy CFL diagnostics.

**Unit Tests**:
- `tests/swe/test_swe_solver.cpp`: Solve known system, verify solution.
- `tests/swe/test_swe_cfl.cpp`: Verify CFL diagnostic calculation without changing `dt`.

---

### 4.4 Task: SWE Source/Sink Terms

**Sub-tasks:**

4.4.1. Implement rainfall (`h += rain_rate * dt`) and evaporation (`h -= evap_rate * dt`).
4.4.2. Implement wind stress (`tau_wind = rho_air * Cw * |W| * W`).
4.4.3. Implement Manning friction (`S_f = g * n^2 * |u| * u / h^(4/3)`).
4.4.4. Implement wet/dry handling (`min_depth = 1e-8`).

**Deliverables**: Source/sink term implementations.

**Unit Tests**: `tests/swe/test_swe_sourcesink.cpp`.

---

### 4.5 Task: SWE Complete Time Stepping

**Sub-tasks:**

4.5.1. Implement `stepForward()` — one complete time step:
1. Apply BCs
2. Assemble RHS (explicit terms)
3. Assemble matrix (implicit terms)
4. Solve (PETSc KSP)
5. Update velocity
6. Wet/dry correction
7. CFL check/diagnostic using legacy behavior; keep the existing Frehg time-loop semantics

4.5.2. Implement `run()` — reproduce Frehg's existing main time loop first, with HDF5 output
at each `dt_out` interval. Do not add SW adaptive stepping until after b1-sw passes.

**Deliverables**: Complete SWE time-stepping loop.

**Unit Tests**:
- `tests/swe/test_swe_timestep.cpp`: One time step, verify mass balance. Multiple steps, verify stability.

---

### 4.6 Validation: Run b1-sw Benchmark

**Procedure**:
1. Run Frehg2 on `benchmarks/b1-sw/b1-sw.yaml`
2. For each output step `t`, compare `depth_*` with legacy reference:
   - Read HDF5 output (Python h5py)
   - Read legacy text output (same Python script as Phase 0.4)
   - Compute L2 error: `||h_new - h_ref||_2 / ||h_ref||_2`
3. Check: L2 error < 1e-6 → **PASS**
4. If FAIL: debug step-by-step
   - Compare matrix assembly (dump intermediate matrices)
   - Compare solver solutions
   - Compare BC application
   - Compare source/sink terms
5. **DO NOT proceed to Phase 5 until this passes.**

---

## 6. Phase 5: Groundwater Module — Predictor-Corrector RE

> Same strategy as Phase 4: port Frehg's Predictor-Corrector RE algorithm exactly, validate with tiny 1D column first, then full b2-gw benchmark.

### 5.0 Task: Tiny Deterministic Test — 1D RE Column

**Sub-tasks:**

5.0.1. Create `tests/re/test_1d_re_column.cpp`:
- 10-layer 1D soil column with uniform VG parameters
- Fixed head BC at top and bottom
- Assemble predictor-corrector matrix
- Dump matrix for predictor step and corrector step
- Compare with legacy Frehg matrix dumps — **must match element-by-element** (tolerance: 1e-12)

**Deliverables**: `tests/re/test_1d_re_column.cpp`.

**Validation**: Matrix element comparison passes.

---

### 5.1 Task: RE Initialization

**Sub-tasks:**

5.1.1. Create `src/re/ReSolver.hpp/cpp`:
- Preserve Frehg's Predictor-Corrector algorithm from `groundwater.c`
- **Do NOT implement** Picard or Newton solvers
- **Do NOT implement** post-allocation step

5.1.2. Implement `initialize()`:
- IC: constant `init_wc` → compute `h` via VG inversion
- BC: fixed head, fixed flux, free outflow
- Read VG parameters from YAML

5.1.3. Implement VG model:
```
Se = 1 / [1 + |alpha*h|^n]^m,  m = 1 - 1/n
theta = theta_r + (theta_s - theta_r) * Se
K = Ks * sqrt(Se) * [1 - (1 - Se^(1/m))^m]^2
```

**Deliverables**: `ReSolver` class with initialization and VG.

**Unit Tests**:
- `tests/re/test_re_init.cpp`: Initialize with constant `wc`, verify VG `theta(h)` and `K(h)`.

**Validation**: Initialize b2-gw benchmark, compare with reference.

---

### 5.2 Task: RE Matrix Assembly

**Sub-tasks:**

5.2.1. Implement `computeKFace()` — use the same arithmetic face averaging as legacy
Frehg's `compute_K_face()` for Part 1 validation.
5.2.2. Implement `assembleMatrixPredictor()` and `assembleMatrixCorrector()`.
5.2.3. Use PETSc `Mat` for matrix assembly.
5.2.4. Index mapping: 3D → 1D PETSc index must match Phase 0 audit.

**Deliverables**: Predictor-Corrector matrix assembly with PETSc.

**Unit Tests**:
- Re-run `tests/re/test_1d_re_column.cpp` — must still pass.
- `tests/re/test_re_matrix.cpp`: Verify coefficients, test arithmetic face averaging.

**Validation**: Dump matrix for b2-gw case, compare with legacy.

---

### 5.3 Task: RE Linear Solver (PETSc)

**Sub-tasks:**

5.3.1. Implement `solveMatrix()` — PETSc KSP.
5.3.2. Implement `updateWaterContent()` — VG model.
5.3.3. Implement `computeDarcyFlux()`.
5.3.4. Implement adaptive time stepping: `dt_GW` from CFL condition.

**Deliverables**: PETSc-based RE solver, adaptive stepping.

**Unit Tests**: `tests/re/test_re_solver.cpp`.

---

### 5.4 Task: RE Source/Sink Terms

**Sub-tasks:**

5.4.1. Infiltration as top BC.
5.4.2. Evapotranspiration as sink term.
5.4.3. Bottom BC: free drainage, fixed head, zero flux.

**Deliverables**: RE source/sink implementations.

**Unit Tests**: `tests/re/test_re_sourcesink.cpp`.

---

### 5.5 Task: RE Complete Time Stepping

**Sub-tasks:**

5.5.1. Implement `stepForward()`:
1. Apply head BCs
2. Predictor: assemble and solve
3. Corrector: update and solve
4. Update water content, Darcy fluxes
5. Adaptive dt check

5.5.2. Implement `run()` — main time loop.

**Deliverables**: Complete RE time-stepping loop.

**Unit Tests**: `tests/re/test_re_timestep.cpp`.

---

### 5.6 Validation: Run b2-gw Benchmark

**Procedure**:
1. Run Frehg2 on `benchmarks/b2-gw/b2-gw.yaml`
2. Compare `head_*`, `moisture_*`, `qz_*` with legacy text reference (HDF5 vs legacy text)
3. L2 error < 1e-6 → **PASS**
4. **DO NOT proceed to Phase 6 until this passes.**

---

## 7. Phase 6: Surface-Subsurface Coupling — Frehg Algorithm

> Phase 6 ports Frehg's **original** coupling algorithm to the new infrastructure. Part 2 Phase 8 will replace it with SERGHEI's asynchronous coupling.

### 6.1 Task: Coupling Mechanism (Frehg Algorithm)

**Sub-tasks:**

6.1.1. Create `src/coupling/Coupling.hpp/cpp`:
- Port Frehg's coupling logic from `solve.c`
- Exchange flux: `q_exchange = K(z_top) * (h_sw - h_gw) / dz`
- Seepage face handling
- Infiltration capacity handling

6.1.2. Implement `computeExchangeFlux()`:
- At each land-surface cell, compute flux between SW and GW
- Same formulation as Frehg

6.1.3. Implement `syncTimeStep()`:
- SW and GW use same `dt` (sync coupling, matching legacy behavior)

**Deliverables**: `Coupling` class with Frehg's original algorithm.

**Unit Tests**:
- `tests/coupling/test_coupling_frehg.cpp`: Known `h_sw`, `h_gw` → verify `q_exchange` matches legacy.
- `tests/coupling/test_coupled_timestep.cpp`: Run coupled case, verify mass balance.

**Validation**:
- Run a coupled case that exercises both SW and GW
- Verify mass balance (SW mass + GW mass + exchange = constant)
- No benchmark comparison at this stage (no coupled reference benchmark in legacy)

---

## 8. Phase 7: Solute Transport Module — Frehg Algorithm

> Phase 7 ports Frehg's **original** solute transport algorithm. Part 2 may extend it.

### 7.1 Task: Solute Transport

**Sub-tasks:**

7.1.1. Create `src/solute/SoluteSolver.hpp/cpp`:
- Port Frehg's solute transport algorithm
- Advection: upwind or TVD scheme
- Diffusion: explicit or implicit
- Use PETSc for implicit diffusion

7.1.2. Implement solute coupling (SW-GW):
- Advective: `Q_exchange * C`
- Diffusive: `D * (C_sw - C_gw) / dz`

**Deliverables**: `SoluteSolver` class.

**Unit Tests**:
- `tests/solute/test_solute_solver.cpp`: Advection of Gaussian pulse, diffusion of step function.
- `tests/solute/test_solute_coupling.cpp`: Verify solute mass balance.

---

### 7.2 Part 1 Completion Checkpoint

Before proceeding to Part 2, verify:
- [ ] `cmake .. && make -j4` — clean build with gcc-15
- [ ] All unit tests pass (`make test`)
- [ ] b1-sw benchmark: L2 error < 1e-6 vs legacy
- [ ] b2-gw benchmark: L2 error < 1e-6 vs legacy
- [ ] Coupling tests pass (mass balance)
- [ ] Solute tests pass
- [ ] HDF5 output files readable and correct

**Part 1 architecture is now verified.** Any bugs from this point are algorithm-related.

---

---

## Part 2: Algorithm Upgrade — New Features from SERGHEI

---

## 9. Phase 8: Asynchronous Coupling — from SERGHEI

> Replace Frehg's synchronous coupling with SERGHEI's asynchronous coupling.

### 8.1 Task: Async Coupling Mechanism

**Sub-tasks:**

8.1.1. Study SERGHEI's coupling implementation (`serghei.h`):
- `async` flag
- Rain → ponding → infiltration logic
- When ponding depth > 0 → infiltration to GW
- GW may discharge to surface (seepage)

8.1.2. Implement `asyncCoupling()` in `Coupling` class:
```
// In SW solver:
// Legacy sign convention: qss < 0 means infiltration (SW -> GW), qss > 0 means seepage (GW -> SW)
h[i] += qss[i] * dt

// In GW solver:
qtop[i] = infiltration_rate  // SW → GW
```

8.1.3. Implement independent time stepping:
- SW and GW can use different `dt`
- GW may take multiple steps per SW step (or vice versa)

8.1.4. Add `async: true` flag to YAML config (default for new simulations).

**Deliverables**: Asynchronous coupling.

**Unit Tests**:
- `tests/coupling/test_async_coupling.cpp`: Rain → ponding → infiltration chain test. Verify mass balance.
- `tests/coupling/test_different_dt.cpp`: SW dt ≠ GW dt, verify correct step count.

**Validation**:
- Run b1-sw and b2-gw benchmarks again — coupling change should NOT affect decoupled benchmarks.
- Run coupled test case — verify physically reasonable exchange flux time series.

---

## 10. Phase 9: Polygon-Based BC & Source/Sink — from SERGHEI

### 9.1 Task: Polygon Class

**Sub-tasks:**

9.1.1. Create `src/bc/Polygon.hpp/cpp`:
- Read polygon from file (text: `x y` per vertex)
- Point-in-polygon test (ray casting)
- `IsInside(x, y)` method

9.1.2. Implement polygon-based BC application:
- Read polygon file referenced in YAML
- For each cell inside polygon → apply specified BC
- Multiple polygons → multiple BC regions

9.1.3. Reference: SERGHEI's `extbc.input` and polygon file handling.

**Deliverables**: `Polygon` class with BC application.

**Unit Tests**:
- `tests/bc/test_polygon.cpp`: Known polygon, test `IsInside()`. Multiple polygons, verify non-overlapping.

**Validation**:
- Run b1-sw with polygon-defined BCs (replacing domain-edge BCs).
- Results should match (if polygon matches domain boundary).

---

### 9.2 Task: Advanced Boundary Conditions

**Sub-tasks:**

9.2.1. Implement all BC types (from SERGHEI + Frehg):
- Fixed water level (Dirichlet, possibly time-varying)
- Fixed flow rate (Neumann, possibly time-varying)
- Free outflow (`dh/dx = 0`)
- Zero-gradient
- Tidal (time series water level)
- Groundwater head
- Root water uptake
- Gravity drainage

9.2.2. Each BC can be applied via polygon OR domain edge.

**Deliverables**: All BC types, polygon-applicable.

**Unit Tests**: `tests/bc/test_advanced_bc.cpp`: Test each BC type.

**Validation**:
- Run benchmarks with new BC system.
- Compare with legacy BC results — should match when same BC applied.

---

### 9.3 Task: Polygon-Based Source/Sink

**Sub-tasks:**

9.3.1. Rainfall, ET, inflow, pumping — apply only inside specified polygon(s).

9.3.2. Multiple source/sink polygons can coexist.

**Deliverables**: Polygon-based source/sink application.

**Unit Tests**: `tests/bc/test_sourcesink.cpp`.

**Validation**: Run b1-sw with polygon-defined rainfall zone — verify only polygon area receives rain.

---

## 11. Phase 10: Non-uniform Soil & Multi-VG — from SERGHEI

### 10.1 Task: Non-uniform Soil Parameters

**Sub-tasks:**

10.1.1. Extend `GwDomain` to support per-cell soil type:
- Read soil map (ASCII raster with soil type IDs)
- `soilID[i]` — soil type index for each cell
- `vgTable[soilID]` — lookup table of VG parameters

10.1.2. Extend YAML config:
```yaml
groundwater:
  soil_map_file: "soil_map.asc"   # ASCII raster with soil IDs
  soil_parameters:
    - id: 0
      name: "sand"
      soil_a: 1.43
      soil_n: 1.56
      wcs: 0.33
      wcr: 0.0
      Ksz: 2.89e-6
    - id: 1
      name: "clay"
      soil_a: 3.5
      soil_n: 1.2
      wcs: 0.45
      wcr: 0.1
      Ksz: 1.0e-7
```

10.1.3. Modify VG model functions to accept `soilID[i]` lookup.

10.1.4. Reference: SERGHEI's `GwDomain.h` — `soilID` and `vgTable` pattern.

**Deliverables**: Non-uniform VG parameter support.

**Unit Tests**:
- `tests/re/test_nonuniform_soil.cpp`: Two soil types, verify different `K` and `theta` at different locations.
- `tests/re/test_soil_map.cpp`: Read soil map ASCII, verify `soilID` assignment.

**Validation**:
- Run b2-gw with uniform soil → same results as Phase 5.
- Run b2-gw with non-uniform soil → results differ from uniform (as expected).
- Verify physically reasonable: sand drains faster than clay.

---

## 12. Phase 11: Advanced Features — Flexible IC, Monitoring Upgrade

### 11.1 Task: Flexible Initial Conditions

**Sub-tasks:**

11.1.1. Support multiple IC types:
- Constant: `init_h: 1.0`
- From file: `init_h_file: ic.h5`
- From variable: `init_wc: 0.3` → compute `h` via VG
- Support multiple variables: `eta`, `h`, `u`, `v` (SW); `h`, `wc`, `psi`, saturation (GW)

11.1.2. IC reader: HDF5, ASCII raster, or constant.

11.1.3. Reference: SERGHEI's flexible IC parser.

**Deliverables**: Flexible IC system.

**Unit Tests**: `tests/core/test_ic.cpp`.

**Validation**: Initialize with file-based IC, compare with constant IC.

---

### 11.2 Task: Monitoring System Upgrade

**Sub-tasks:**

11.2.1. Polygon-based monitoring:
- User specifies polygon for monitoring region
- Model outputs all variables at all cells inside polygon
- Output: `monitor_polygon_001.h5`

11.2.2. Keep legacy point monitoring for backward compatibility.

11.2.3. Mass balance monitoring: total mass, boundary fluxes, source/sink.

**Deliverables**: Upgraded monitoring system.

**Unit Tests**: `tests/core/test_monitor.cpp`.

**Validation**: Run benchmark with polygon monitoring, verify output.

---

### 11.3 Task: Performance Optimization

**Sub-tasks:**

11.3.1. Kokkos tuning:
- Appropriate execution space for each kernel
- Minimize host-device copies

11.3.2. PETSc tuning:
- Optimal solver/preconditioner per equation type
- Command-line options for runtime tuning

11.3.3. Memory optimization: use `float` where possible.

**Deliverables**: Optimized code.

**Validation**: Performance benchmark — Frehg2 should be faster than Frehg.

---

## 13. Phase 12: Validation, Documentation & Release

### 12.1 Task: Complete Benchmark Validation

**Sub-tasks:**

12.1.1. Run all benchmarks:
- `b1-sw`: SW only — must match legacy
- `b2-gw`: GW only — must match legacy
- Coupled case: verify mass balance

12.1.2. Run new feature tests:
- Polygon BC test (SERGHEI-style)
- Non-uniform soil test
- Flexible IC test
- Async coupling test

**Deliverables**: Validation report.

---

### 12.2 Task: Documentation

**Sub-tasks:**

12.2.1. User manual: installation, YAML config reference, input formats, running benchmarks.
12.2.2. Developer manual: code architecture, adding modules, testing.
12.2.3. Theory manual: governing equations, numerical methods.

---

### 12.3 Task: Release

**Sub-tasks:**

12.3.1. Tag version `v2.0.0`.
12.3.2. Create release notes.
12.3.3. Archive legacy code.

---

## Summary

### Part 1: Architecture (Phases 0-7)

| Phase | Task | Key Validation |
|-------|------|---------------|
| 0 | Legacy audit | All BC codes, indices, variables, formats documented |
| 1 | Build system | CMake finds Kokkos/MPI/PETSc, gcc-15 compiles |
| 2 | Data structures | Grid, Domain, State with Kokkos Views |
| 3 | I/O layer | YAML config parses, HDF5 writes/reads, ASCII raster reads |
| 4 | SWE solver | 1D matrix test → **b1-sw benchmark matches legacy** |
| 5 | RE solver | 1D column test → **b2-gw benchmark matches legacy** |
| 6 | Coupling (Frehg) | Mass balance verification |
| 7 | Solute (Frehg) | Analytical solution comparison |

### Part 2: Algorithm (Phases 8-12)

| Phase | Task | Key Validation |
|-------|------|---------------|
| 8 | Async coupling | Rain→ponding→infiltration chain, b1-sw/b2-gw still pass |
| 9 | Polygon BC/source | Polygon BC matches domain-edge BC for same inputs |
| 10 | Non-uniform soil | VG parameters per soil type, physically reasonable drainage |
| 11 | Flexible IC, monitoring | Multiple IC types, polygon monitoring output |
| 12 | Validation, release | All benchmarks pass, docs complete, v2.0.0 tagged |

> **CRITICAL RULES**:
> 1. Part 1 preserves Frehg's algorithms exactly — validated against legacy output.
> 2. Part 2 adds new features incrementally — re-run benchmarks after each change.
> 3. Never mix architecture changes with algorithm changes in the same step.
> 4. Use gcc-15 only; dependencies at `/Users/zhili/Codes/local/` (do NOT re-install).
> 5. Each tiny deterministic test (1D SWE/RE) must pass before running full benchmarks.

---

## 14. Appendix: File Format Specifications

### A.1 YAML Configuration File Format (Frozen after Phase 0)

```yaml
simulation:
  title: "Benchmark 1: Surface Water Only"
  author: "Your Name"
  code_version: "2.0.0"

domain:
  nx: 1
  ny: 10
  nz: 1
  dx: 80.0
  dy: 80.0
  dz: 0.1
  dz_incre: 1.0
  botZ: -3.0
  mpi_nx: 1
  mpi_ny: 1

time:
  dt: 5.0
  Tend: 18000.0
  dt_out: 1800.0
  cfl: 0.9
  Co_max: 2.0

surface_water:
  enable: true
  solver: "semi_implicit"
  init_eta: -2.0
  bc_type: [0, 0, 0, 0]   # Legacy bctype_SW order from Phase 0 audit
  rain_file: "rain.csv"
  manning: 0.019
  min_depth: 1.0e-8
  grav: 9.81

groundwater:
  enable: true
  solver: "predictor_corrector"
  init_wc: 0.033
  use_vg: true
  soil_a: 1.43
  soil_n: 1.56
  wcs: 0.33
  wcr: 0.0
  Ksz: 2.89e-6
  bc_type: [0, 0, 0, 0, 0, 0]   # Legacy bctype_GW: x+, x-, y+, y-, bottom(z+), top(z-)
  dt_min: 1.0e-4
  dt_max: 2.0

solute:
  enable: false

output:
  format: "hdf5"
  filename: "output.h5"
  variables:
    - "water_depth"
    - "velocity_x"
    - "velocity_y"

monitor:
  points:
    - name: "point_001"
      locX: 5
      locY: 5
      variables: ["water_depth", "velocity_x"]
```

### A.2 Polygon File Format

```
# format: x y (one vertex per line, close polygon by repeating first point)
0.0 0.0
100.0 0.0
100.0 100.0
0.0 100.0
0.0 0.0
```

### A.3 Time Series File Format (CSV)

```csv
time,value
0.0,0.0
3600.0,0.001
7200.0,0.002
```

### A.4 ASCII Raster Format (ESRI)

```
ncols 100
nrows 100
xllcorner 0.0
yllcorner 0.0
cellsize 80.0
NODATA_value -9999
100.0 101.0 102.0 ...  (data rows, bottom to top)
```
