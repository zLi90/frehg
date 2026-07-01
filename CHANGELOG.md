# Changelog

All notable changes to Frehg2 are documented in this file. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-06-30

First production release. Frehg2 is a complete C++20 / Kokkos / MPI / CMake / PETSc rewrite of the
legacy serial C / MPI / LASPack `frehg` code, delivered phase-by-phase per `INTEGRATED_PLAN.md`
(archived under `docs/planning/`) with the legacy algorithm reproduced to gate tolerance.

### Added

- **Build & core infrastructure (P1–P2).** Reproducible CMake build (`g++-15`, C++20, strict
  warnings-as-errors); core types (`real`, Kokkos `View` aliases); halo-padded flat `Grid`;
  model-owned `DomainDecomposition` (5-/7-point) and hand-written `MpiComm` halo exchange — **no
  PETSc `DMDA`**.
- **Backend-agnostic linear solver (P2/P2B).** `LinearSolver`/`SparseSystem` seam; PETSc KSP
  backend (`MATMPIAIJ`, COO assembly) is the only PETSc-linked library; physics never sees a
  solver-library type. Solver bake-off harness.
- **I/O layer (P3).** YAML `Config`; `OutputWriter` seam with HDF5 (C API) default, three
  `io_mode`s, XDMF sidecar, embedded provenance (git SHA, config SHA-256, units), and atomic
  checkpoint/restart.
- **Shallow-water solver (P4).** Legacy-exact `b1-sw` port; well-balanced lake-at-rest (`b0-lake`);
  serial + MPI rank-equivalence gates.
- **Richards solver (P5).** Legacy-exact predictor-corrector RE with Van Genuchten constitutive
  relations and 3-criterion adaptive time stepping; element-wise legacy parity (`b2-gw`).
- **SW↔GW coupling (P6, P11).** Conservative synchronous (Frehg-style) coupling, plus an opt-in
  double-buffered async pipeline (`coupling.mode: async`) that is bit-identical to the sequential
  path.
- **Orchestrator (P7).** The single general-purpose production driver — the only path from `main()`
  to physics; coupled time-stepping, full-state checkpoint/restart, `simulation_summary.txt`.
- **Solute transport (P8, P16).** Advection (upwind / MUSCL), implicit diffusion (through the seam),
  source/sink mixing, first-order decay; integrated into the Orchestrator under `solute.enabled`.
- **On-node Kokkos + GPU backend (P9, P10).** Single numerical-free Kokkos pass over every per-cell
  loop (bit-identical CPU/OpenMP); GPU-capable Kokkos→PETSc COO assembly (`MatSetValuesCOO`) and
  guarded GPU async solver groups behind compile guards.
- **Polygon BC & sources (P12).** Polygon-based boundary conditions and source/sink regions
  (discharge, depth, critical-weir, inflow, rainfall, extraction wells).
- **Non-uniform soil (P13).** Per-cell soil class map (CSV / ESRI raster) with per-class VG
  parameters.
- **Flexible initial conditions (P14).** Constant / raster / polygon / formula / restart ICs.
- **Monitoring (P15).** Point probes, line-flux integration, CSV monitor output.
- **YAML schema v2 + migration (P17).** Documented `docs/yaml_schema_v2.md` (generated), v1→v2
  migration tool, and an Orchestrator schema-version gate.
- **Benchmarks b3–b6 (P18).** Review-tier ports with per-benchmark docs and capped-step integration
  tests.
- **Unified validation (P19).** `tools/run_validation.py` three-tier (PASS/REVIEW/FAIL) harness +
  real water mass-balance budget in the run summary; CI workflow.
- **Performance instrumentation (P21).** `Frehg2::perf` (Kokkos timers + profiling regions),
  `perf_*` breakdown in `simulation_summary.txt`, `tools/perf_report.py`, and three staging/setup
  optimizations.
- **Release docs (P22).** `README.md`, `docs/architecture.md`, `docs/developer_guide.md`,
  `docs/coding_standards.md` (+ `.clang-format`), `docs/roadmap.md`, AMR/MG research notes,
  `CHANGELOG.md`, `CITATION.cff`, `LICENSE` (BSD 3-Clause).

### Changed

- **Legacy `frehg` deprecated (P20).** The legacy C/MPI/LASPack code is archived under `legacy/`,
  never compiled by default, and never linked or `#include`d by any Frehg2 target. Opt-in via
  `FREHG2_USE_LEGACY=ON` (forwards to its own historical makefile). See
  `docs/migration_from_legacy.md`.

### Validation status

- `b1-sw` **PASS** (legacy depth L2 = 0 vs the committed reference); `b0-lake` **PASS**
  (`max|eta-1| ~ 2e-15`).
- `b2-gw`, `b3-kirkland`, `b4-govindaraju`, `b5-vcatchment`, `b6-kuan` **REVIEW** (registry-`review`
  gate; stable, finite, sound mass balance/trends). Overall **REVIEW** — as expected, five of six
  are review-tier by design.

### Known limitations / deferred

- **GPU execution** is validated later on a Linux/NVIDIA host (`-DFREHG2_ENABLE_CUDA=ON`); GPU tests
  are `gpu`-labeled and DISABLED on macOS. See `docs/gpu_validation_policy.md`.
- **MPI strong-scaling** (≥ 0.7 efficiency at 16 ranks) and **GPU speedup** (≥ 5× at 1M cells) are
  deferred to the HPC/CUDA host; instrumentation is in place.
- **Richards** is vertical-only with a per-column soil map (no lateral/z-layered RE); the full 2-D
  `b3-kirkland` ships as a 1-D vertical surrogate.
- **Signed `v1.0.0` git tag + GitHub release** are pending repository `git init` (the repo is not
  yet under version control); build provenance auto-records the git SHA once it is.

[1.0.0]: https://github.com/frehg2/frehg2/releases/tag/v1.0.0
