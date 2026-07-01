# Frehg2 Architecture

Frehg2 is a surface-water (SWE) — groundwater (Richards) — solute-transport coupled model. It is
built as a set of small static libraries with a strict dependency order, a single
backend-agnostic linear-solver seam, and **exactly one production code path**: `Orchestrator::run()`.

This document is the high-level map. For build/test mechanics see
[`developer_guide.md`](developer_guide.md); for the YAML configuration see
[`yaml_schema_v2.md`](yaml_schema_v2.md).

---

## 1. Design principles (recap)

These govern every module (full text in `INTEGRATED_PLAN.md`, "Design Principles DP1–DP9"):

- **DP1 — One Orchestrator, one production path.** Benchmarks, validation, and user runs all go
  through `Orchestrator::run()`. There is no `main_production.cpp`, no `runProductionCoupled()`.
- **DP2 — Backend-agnostic linear algebra.** All global solves go through the
  `LinearSolver` / `SparseSystem` interface. Physics code never sees `Mat`, `Vec`, `KSP`, `DM`,
  or any solver-library type. The model owns its parallel layout (`DomainDecomposition`) and halo
  exchange (`MpiComm`) — there is **no PETSc `DMDA`**.
- **DP3 — Kokkos for local compute, interface for global assembly; numerical-first, Kokkos/GPU-last.**
  Per-cell kernels are Kokkos `parallel_for`; assembly emits COO triplets to the seam. The full
  numerical model was built and validated as plain serial/MPI loops first (through P16); the Kokkos
  pass (P9) and GPU pass (P10) converted the entire kernel surface at once with **no new numerical
  behavior**.
- **DP4 — Variable convention locked.** Surface `eta, u, v`; subsurface `pressure_head`; solute
  `conc`. The surface matrix is `A · eta_new = b`. (Not the SERGHEI `h, hu, hv`.)
- **DP5 — Arithmetic-mean K-face**, `0.5·(Kp + Km)`.

---

## 2. Library / module map

Each directory under `src/` builds one static library with a `Frehg2::<name>` alias. Public headers
live under `include/frehg2/<module>/`.

| Library | Alias | Responsibility | Links (direct) |
|---------|-------|----------------|----------------|
| `frehg2_core` | `Frehg2::core` | `Grid` (halo-padded flat storage), `State`/`GwState`, `MpiComm` halo exchange, `ParallelFor` wrappers, `define.hpp`/`types.hpp` | Kokkos, MPI |
| `frehg2_decomp` | `Frehg2::decomp` | `DomainDecomposition` (`Decomp2D` 5-pt / `Decomp3D` 7-pt), global row numbering. **No PETSc.** | core |
| `frehg2_petsc_backend` | `Frehg2::petsc_backend` | The **only** PETSc-linked library: `PetscLinearSolver`, `PetscSparseSystem`, Kokkos→PETSc COO bridge (host Option A / device Option B), `makeLinearSolver` factory, GPU subcomm split | decomp, PETSc |
| `frehg2_io` | `Frehg2::io` | `Config` (YAML), `OutputWriter` seam + HDF5 writer/reader, XDMF sidecar, checkpoints, provenance, SHA-256, schema migration | core, yaml-cpp, HDF5 |
| `frehg2_swe` | `Frehg2::swe` | Shallow-water solver (`SweSolver`), legacy-exact `b1-sw` path | core (+ seam) |
| `frehg2_re` | `Frehg2::re` | Richards solver (`ReSolver`), Van Genuchten constitutive, predictor-corrector + adaptive dt | core, soil (+ seam) |
| `frehg2_soil` | `Frehg2::soil` | `SoilMap` (per-cell class table + raster/CSV), `SoilParams` | core, io |
| `frehg2_solute` | `Frehg2::solute` | `SoluteStepper`: advection (upwind/MUSCL), implicit diffusion (seam), source/sink, decay | core, io (+ seam) |
| `frehg2_coupling` | `Frehg2::coupling` | Synchronous SW↔GW exchange (`Coupling`), conservative apply | core, swe, re |
| `frehg2_bc` | `Frehg2::bc` | Polygon BC + source/sink regions (ray-casting, `PolygonIndex`) | core, io |
| `frehg2_ic` | `Frehg2::ic` | Flexible initial conditions: constant/raster/polygon/formula/restart | core, io, swe, re, bc |
| `frehg2_monitoring` | `Frehg2::monitoring` | Point probes, line-flux integration, CSV monitor output | core, io |
| `frehg2_perf` | `Frehg2::perf` | `ScopedTimer`, `PerfRecorder`, `Counters` (Kokkos timers + profiling regions) | Kokkos, MPI |
| `frehg2_orchestrator` | `Frehg2::orchestrator` | The production driver: `Orchestrator`, `AsyncPipeline`, `DoubleBuffer` | coupling, bc, ic, monitoring, solute, io, perf, petsc_backend |

The `frehg2` executable (`src/driver/main.cpp`) links only `Frehg2::orchestrator`, which pulls the
rest transitively.

---

## 3. Dependency graph

```
                 Kokkos + MPI
                      │
                   core ──────────────► perf
                    │  │
        ┌───────────┘  └───────────┐
      decomp                      io (yaml-cpp, HDF5)
        │                          │
  petsc_backend ◄──────────┐   ┌───┴─────────┬──────────┐
   (ONLY PETSc lib)        │  soil          bc       monitoring
        ▲                  │   │             │
        │ (LinearSolver    │   re ◄──────────┘
        │  seam, no PETSc  │   │
        │  in physics)     │  swe
        │                  │   │
        │              coupling (swe + re)
        │                  │
        │                 ic (swe, re, bc)
        │                  │
        └────────── orchestrator ──────────► frehg2 (main.cpp)
```

Key invariants (enforced by CMake seam-check targets — see §5):

- **PETSc appears only in `frehg2_petsc_backend`.** `decomp` is provably PETSc-free
  (`test_decomp_no_petsc`).
- Physics libs (`swe`, `re`, `solute`, `coupling`) talk only to the
  `LinearSolver`/`SparseSystem` interface headers; they never include `petsc.h`.

---

## 4. The production path (DP1)

```
main()                                  src/driver/main.cpp
  └─ Config::fromFile(yaml)             src/io/Config.cpp
  └─ Orchestrator::initialize(config)   src/driver/Orchestrator.cpp
        ├─ validateSchemaV2()           (refuse schema_version != "2.0")
        ├─ build Grid / MpiComm / DomainDecomposition
        ├─ makeLinearSolver()           backend-agnostic factory (PETSc inside)
        ├─ build enabled solvers        SweSolver / ReSolver per modules.*
        ├─ build soil / IC / BC / monitors / solute / output
        └─ setupInitialConditions()
  └─ Orchestrator::run()
        └─ loop: step() → output/checkpoint/monitors at cadence
```

`Orchestrator::step()` advances one coupling window with synchronous, fully-advanced semantics:

```
step():
  swAdvanceTo(t_new)              SWE solve + rain/evap + polygon surface BC/sources
  while gw_time < t_new:          GW adaptive substeps catch up to SW time
      gwCatchUpTo(t_new)          RE predictor-corrector + polygon wells
  applyCouplingExchange()         conservative SW↔GW top-face exchange
  applySolute()                   operator-split transport pass (if solute.enabled)
```

### Async coupling (P11)

`coupling.mode: async` (single rank only) runs the **identical** Gauss–Seidel sequence but overlaps
the groundwater catch-up of the *previous* window (on a worker thread) with the SW advance of the
current window, via `src/driver/{DoubleBuffer.hpp, AsyncPipeline.cpp}`:

```
window k:   SW advance (caller thread)        ║  GW catch-up of window k-1 (worker thread)
            ─────────────────────────────────╫──────────────────────────────────────────
            join ──► applyExchange(window k)  ║
```

Because the two threads touch disjoint `SweFields`/`GwFields` and both depend only on the
already-applied exchange, the async result is **bit-identical** to the sequential path (gated at
`< 1e-12`). PETSc KSP solves are serialized (shared communicator is not thread-safe); multi-rank
async falls back to synchronous. True concurrent GPU solves (PetscSubcomm + CUDA streams) are
designed in [`research_notes/async_gpu.md`](research_notes/async_gpu.md) and deferred to the GPU
validation host.

---

## 5. The linear-solver seam (DP2/DP3)

```
physics kernel (Kokkos)                  src/{swe,re,solute}
   │  per-cell stencil coefficients
   ▼
SparseSystem::addRow / addCOO            include/frehg2/linear/SparseSystem.hpp
   │  COO triplets (global_row, global_col, value)
   ▼
PetscSparseSystem  (backend only)        src/linear/backends/
   │  CPU: MatSetValues on HostSpace (Option A)
   │  GPU: MatSetPreallocationCOO + MatSetValuesCOO (Option B, FREHG2_GPU_ASSEMBLY)
   ▼
LinearSolver::solve(A, b, x)             PetscLinearSolver (KSP; cg/gmres + pc)
```

Global row/column indices come from `DomainDecomposition::globalRow` / `::stencilColumns`. Matrices
are `MATMPIAIJ` (or `MATAIJKOKKOS` on GPU) preallocated from the model decomposition — never a
`DMDA`, never `MatCreateSeqAIJ` (serial = one-rank MPIAIJ).

**Static enforcement** (CMake custom targets, run as build dependencies of `frehg2`):

- `check_no_seqaij` — no `MatCreateSeqAIJ` / `PETSC_COMM_SELF` / `DMDA*` in `src/`.
- `check_solver_seam` — no `MatSetValues` outside the bridge; no `Mat`/`Vec`/`KSP`/`DM`/`Tpetra` in
  `src/{swe,re,solute,coupling}`.
- `check_gpu_coo_assembly` — the device assembly TU uses `MatSetValuesCOO` only.

---

## 6. Data layout & parallelism

- **Storage:** flat, **halo-padded** `Kokkos::View`s. One ghost on each horizontal side; the
  subsurface is halo-padded in z too. `Grid::getIndex(i,j,k) = (i+1) + (j+1)·(nx+2) + (k+1)·(nx+2)·(ny+2)`.
- **Decomposition:** x/y only; z is on-rank (no MPI message in z). Global numbering is contiguous
  within a rank, rank-major — matching PETSc's Vec/Mat layout exactly.
- **Halo exchange:** hand-written `MPI_Sendrecv` in `MpiComm` (`haloExchange2D`/`haloExchange3D`).
  Any kernel reading a neighbor cell must exchange first.
- **On-node:** every per-cell loop goes through `include/frehg2/core/ParallelFor.hpp`
  (`parallelForRange`/`Surface`/`Volume`), `LoopExec = DefaultHostExecutionSpace` (OpenMP here,
  device on a CUDA build). Reductions use `Kokkos::Max`/`Min` for order-independent (bit-identical)
  results. Loops kept deliberately sequential (Gauss–Seidel sweeps, scatter limiters, COO assembly)
  are catalogued in [`local_loops_audit.md`](local_loops_audit.md).

---

## 7. I/O & provenance

- `OutputWriter` seam; HDF5 is the default backend. Fields are written as owned cells (no halo) +
  a global-index map so any decomposition reassembles to the same on-disk global field.
- `output.io_mode`: `serial_gather` (default, gzip), `parallel_collective` (MPI-IO, uncompressed),
  or `file_per_rank`.
- Every file embeds provenance in `/simulation/metadata`: `git_sha`, `git_dirty`, `config_sha256`,
  build type/compiler, Kokkos backend, solver backend, MPI ranks, io_mode; each field carries a
  `units` attribute. An XDMF `.xmf` sidecar is emitted for ParaView/VisIt.
- Checkpoints are atomic (temp → `H5Fflush` → `fsync` → `rename(2)`), one file per set, pruned to
  `time.max_checkpoints`. Restart restores the full halo-padded solver state (serial-only).

---

## 8. Performance instrumentation (P21)

`Frehg2::perf` wraps named regions with `Kokkos::Timer` + optional `Kokkos::Profiling` push/pop, so
the same regions show up in nvprof/VTune on the GPU host. `simulation_summary.txt` carries
`perf_*_seconds` per region, per-step averages, and counters (`perf_cells_touched`,
`perf_ksp_iterations`, `perf_bytes_staged`). `tools/perf_report.py` runs scaling sweeps; hot-spot
analysis is in [`perf/p21_hotspots.md`](perf/p21_hotspots.md).
