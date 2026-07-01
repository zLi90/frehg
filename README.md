# Frehg2

Frehg2 is a production-grade, general-purpose, GPU-capable, C++20 / Kokkos / MPI /
CMake / PETSc surface-water — groundwater — solute-transport coupled numerical
model. It is the modern successor to the legacy serial C / MPI / LASPack `frehg`
code, rebuilt without losing algorithmic fidelity.

The complete execution plan lives in [`INTEGRATED_PLAN.md`](INTEGRATED_PLAN.md);
the engineering rules every change must satisfy live in
[`.cursorrules`](.cursorrules). Read both before contributing.

> **The legacy `frehg` code is deprecated (Phase 20).** It is archived under
> [`legacy/`](legacy/README.md) for reproducibility/citation and is **not compiled by default**.
> Frehg2 (this repository) is the only supported, default build. To translate a legacy `input`
> file to the Frehg2 v2 YAML schema, see [`docs/migration_from_legacy.md`](docs/migration_from_legacy.md).

## Repository layout

| Path                  | Purpose                                                              |
|-----------------------|---------------------------------------------------------------------|
| `INTEGRATED_PLAN.md`  | Authoritative phase-by-phase upgrade plan                            |
| `.cursorrules`        | Mandatory engineering rules / anti-shortcut policy                   |
| `CMakeLists.txt`      | Root build (Frehg2 only; legacy is opt-in via `FREHG2_USE_LEGACY`)   |
| `src/`                | Frehg2 implementation (driver, solvers, I/O, coupling, …)           |
| `include/frehg2/`     | Public headers (`include/frehg2/**/*.hpp`)                           |
| `tests/`              | CTest suite (unit, integration, validation)                         |
| `tools/`              | Tooling: validation harness, schema doc gen, bake-off, migration    |
| `scripts/`            | Legacy-comparison + legacy→YAML conversion scripts                   |
| `benchmarks/`         | Frehg2 benchmark configs `b0-lake … b6-kuan`                         |
| `docs/`               | Schema, validation, benchmarks, migration, legacy audit             |
| `legacy/`             | **Archived** reference material (read-only); see `legacy/README.md` |

## Toolchain (this development machine: macOS / Darwin)

- Compiler: `gcc-15` / `g++-15` (Homebrew)
- MPI: local MPICH (`/Users/zhili/Codes/local/bin/{mpicxx,mpicc,mpiexec}`), pinned to the build
  PETSc was compiled against
- CMake ≥ 3.20
- Dependency prefix: `/Users/zhili/Codes/local/` (PETSc, Kokkos, HDF5, yaml-cpp)
- GPU: **not** available on macOS (CUDA-only); GPU code is written to be GPU-capable but its
  execution is validated later on a Linux/NVIDIA host (`-DFREHG2_ENABLE_CUDA=ON`).

## Build & test

```bash
cmake -S . -B build -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/frehg2 benchmarks/b1-sw/b1-sw.yaml
```

Multi-rank runs use the local MPICH `mpiexec`:

```bash
/Users/zhili/Codes/local/bin/mpiexec -n 2 ./build/frehg2 benchmarks/b1-sw/b1-sw.yaml
```

## Validation

The unified three-tier validation harness is the canonical "is Frehg2 working?" check:

```bash
python tools/run_validation.py --bin build/frehg2 --out build/validation   # all 6 benchmarks
python tools/run_b0_lake.py    --bin build/frehg2                          # well-balanced check
```

See [`docs/validation.md`](docs/validation.md) for the PASS/REVIEW/FAIL tiers and reports.

## Documentation

| Document | Purpose |
|----------|---------|
| [`docs/user_manual.md`](docs/user_manual.md) | **New users start here** — build, configure, run, and post-process; full configuration field reference and sample inputs |
| [`docs/architecture.md`](docs/architecture.md) | Module map, library dependency graph, the single Orchestrator production path, async-coupling pipeline |
| [`docs/developer_guide.md`](docs/developer_guide.md) | Build options, Kokkos/PETSc configuration, test layout, how to add a benchmark or a module |
| [`docs/coding_standards.md`](docs/coding_standards.md) | Naming, formatting (`.clang-format`), CQ1–CQ6 enforcement |
| [`docs/yaml_schema_v2.md`](docs/yaml_schema_v2.md) | Authoritative v2 YAML configuration reference (generated) |
| [`docs/benchmarks/`](docs/benchmarks/) | Per-benchmark validation docs (`b3`–`b6`) and three-tier thresholds |
| [`docs/validation.md`](docs/validation.md) | Unified three-tier validation harness |
| [`docs/perf/p21_hotspots.md`](docs/perf/p21_hotspots.md) | Performance instrumentation + hot-spot analysis |
| [`docs/roadmap.md`](docs/roadmap.md) | Post-1.0 plan (1.1 / 1.2 / 2.0) |
| [`docs/research_notes/`](docs/research_notes/) | AMR, multigrid, multi-GPU, async-GPU design notes |
| [`docs/planning/`](docs/planning/) | Archived `INTEGRATED_PLAN.md` (the executable upgrade plan) |

## Citation

If you use Frehg2 in academic work, please cite it (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@software{frehg2_2026,
  title        = {Frehg2: a production-grade surface-water--groundwater--solute
                  coupled numerical model},
  author       = {{The Frehg2 Development Team}},
  year         = {2026},
  version      = {1.0.0},
  url          = {https://github.com/frehg2/frehg2},
  note         = {C++20 / Kokkos / MPI / PETSc successor to the legacy Frehg code}
}
```

## License

Frehg2 is released under the BSD 3-Clause License — see [`LICENSE`](LICENSE).

## Legacy code (archival, opt-in)

The legacy model is **not** part of the Frehg2 build. To build it for archival/reproducibility:

```bash
cmake -S . -B build-legacy -DFREHG2_USE_LEGACY=ON && cmake --build build-legacy --target legacy_frehg
# or directly with its own historical build system:
make -C legacy/frehg
```

See [`legacy/frehg/DEPRECATED.md`](legacy/frehg/DEPRECATED.md).
