# Frehg2 Developer Guide

How to build, test, and extend Frehg2. For the high-level design see
[`architecture.md`](architecture.md); for style rules see [`coding_standards.md`](coding_standards.md).

---

## 1. Toolchain

This is the development machine (macOS / Darwin, Apple Silicon):

- **Compiler:** `gcc-15` / `g++-15` (Homebrew), C++20 (`CMAKE_CXX_STANDARD 20`, no extensions).
- **MPI:** local MPICH at `/Users/zhili/Codes/local/bin/{mpicxx,mpicc,mpiexec}`, pinned to the build
  PETSc was compiled against. **Use this `mpiexec` for multi-rank runs** — Homebrew's `mpiexec`
  launches the MPICH-linked binary as size-1 singletons.
- **Dependencies** (all under `/Users/zhili/Codes/local/`): PETSc 3.25.1 (pkg-config), Kokkos 5.1.1
  (`OpenMP;Serial`, `Kokkos_CXX_STANDARD=20`), HDF5 (C API, parallel-enabled), yaml-cpp.
- **GPU:** not available on macOS (CUDA only). GPU code is written but its execution is validated
  later on a Linux/NVIDIA host. See [`gpu_validation_policy.md`](gpu_validation_policy.md).

---

## 2. Build

```bash
cmake -S . -B build \
  -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

### CMake options

| Option | Default | Effect |
|--------|:-------:|--------|
| `BUILD_TESTING` | `ON` | Build and register the CTest suite |
| `FREHG2_STRICT_WARNINGS` | `ON` | `-Wall -Wextra -Werror -Wno-unused-parameter` (warnings are errors) |
| `FREHG2_ENABLE_CUDA` | `OFF` | Enable the CUDA GPU backend (Linux/NVIDIA only; **FATAL on macOS**); un-disables `gpu` tests |
| `FREHG2_USE_GPU_MAT` | `OFF` | Default to `MATAIJKOKKOS`/`VECKOKKOS` PETSc types (override at runtime with `-mat_type`) |
| `FREHG2_USE_LEGACY` | `OFF` | Build the archived legacy `frehg` via its own makefile (opt-in; never linked by Frehg2) |

### Runtime solver options (PETSc)

The linear backend is selected through PETSc options, not code:

```bash
./build/frehg2 config.yaml -ksp_type cg -pc_type jacobi -ksp_rtol 1e-12
```

Defaults: `cg` for SPD systems, `gmres` for non-symmetric. Locally there is no direct LU for
`mpiaij` (no MUMPS/SuperLU), so deterministic serial bring-up uses `cg + jacobi` at tight `rtol`;
`-pc_type lu` works unchanged on a MUMPS-enabled (Linux) build.

### Multi-rank

```bash
/Users/zhili/Codes/local/bin/mpiexec -n 2 ./build/frehg2 benchmarks/b1-sw/b1-sw.yaml
```

Serial is one rank (`mpirun -n 1`) — there is no separate sequential code path.

---

## 3. Test layout

Tests mirror `src/` (`tests/<module>/test_*.cpp`) plus `tests/integration/` and `tests/validation/`.
The harness is a vendored minimal Catch2-style header `tests/frehg2_test.hpp` (no Catch2/GTest
installed locally; keep test sources Catch2-API-compatible).

### CTest labels

| Label | Meaning |
|-------|---------|
| `cpu` | Runs on CPU/OpenMP (default) |
| `mpi` | Multi-rank (`mpiexec -n 2/4`) |
| `integration` | Multi-module / Orchestrator path |
| `regression` | Numerical regression vs a reference |
| `validation` | The unified harness checks (P19) |
| `perf` | Perf recorder + report smoke (P21) |
| `gpu` | Requires CUDA + Kokkos-aware PETSc; **DISABLED locally**, deferred to Linux/NVIDIA |
| `benchmark` | Capped-step benchmark integration tests (b3–b6) |

```bash
ctest --test-dir build -R test_swe_b1 --output-on-failure   # one test
ctest --test-dir build -L mpi --output-on-failure           # by label
```

### Seam / anti-pattern checks (build-time)

```bash
cmake --build build --target check_no_seqaij check_solver_seam check_gpu_coo_assembly
```

These also run automatically as build dependencies of the `frehg2` target. They enforce DP2/DP3:
no `MatCreateSeqAIJ`/`PETSC_COMM_SELF`/`DMDA` in `src/`, no solver-library types in physics, and
COO-only GPU assembly.

### Memory safety

AddressSanitizer is unavailable on this toolchain (Homebrew gcc on Darwin has no `libasan`;
Valgrind is unavailable on Apple Silicon). Memory safety is **structural**: `src/` uses only
`Kokkos::View` + PETSc RAII destructors — no raw `new`/`delete`/`malloc`. On a Linux validation
build, run:

```bash
cmake -S . -B build-asan -DCMAKE_CXX_FLAGS="-fsanitize=address" && cmake --build build-asan
ctest --test-dir build-asan --output-on-failure
```

---

## 4. Validation & performance

```bash
python tools/run_validation.py --bin build/frehg2 --out build/validation   # all 6 + b0-lake
python tools/run_validation.py --quick                                     # fast CI smoke
python tools/run_b0_lake.py    --bin build/frehg2                          # well-balanced 1e-12
python tools/perf_report.py    --bin build/frehg2 --out build/perf         # scaling sweep
python tools/perf_report.py    --quick                                     # CI smoke
```

See [`validation.md`](validation.md) and [`perf/p21_hotspots.md`](perf/p21_hotspots.md).

---

## 5. How to add a new benchmark

1. Create `benchmarks/<id>/<id>.yaml` (v2 schema; see [`yaml_schema_v2.md`](yaml_schema_v2.md)) and
   copy any forcing artifacts (DEM `.asc`, `bath`, `rain` series in m/s) beside it.
2. Register it in [`../benchmarks/reference_registry.yaml`](../benchmarks/reference_registry.yaml)
   with its gate (`pass`/`review`) and tolerance. The registry is **authoritative** — it overrides
   plan prose if they disagree.
3. Add an integration test `tests/integration/test_<id>.cpp` (use `tests/integration/benchmark_util.hpp`).
   Cap long runs with `time.max_steps`. Assert *real physics* (finite/non-negative state, mass
   conservation, bounded `wc ∈ [θr, θs]`) — never "no crash".
4. Add per-benchmark thresholds to [`../tools/validation_thresholds.yaml`](../tools/validation_thresholds.yaml)
   and a doc under [`benchmarks/`](benchmarks/) with the three-tier (PASS/REVIEW/FAIL) bands.
5. Confirm it appears in `python tools/run_benchmark.py --list` and passes the harness.

---

## 6. How to add a new module

1. Create `src/<module>/` with a `CMakeLists.txt` defining `frehg2_<module>` + a `Frehg2::<module>`
   alias; add `add_subdirectory(<module>)` in `src/CMakeLists.txt`. Apply `apply_strict_warnings()`.
2. Public API goes in `include/frehg2/<module>/*.hpp` (header guards, **not** `#pragma once` — CQ2);
   implementation in `src/<module>/*.cpp` (every `.hpp` has a `.cpp` — CQ3).
3. **Linear algebra goes through the seam only** — include
   `include/frehg2/linear/{SparseSystem,LinearSolver,SolverConfig}.hpp`, never `petsc.h`. PETSc may
   appear only under `src/linear/backends/`. The seam checks will fail the build otherwise.
4. **Reachability (DP1):** wire the module into `Orchestrator` so there is a path from `main()` to
   every public method. Test-only wiring is a bug.
5. **Numerical-first (DP3):** implement and validate as plain serial/MPI loops; on-node Kokkos
   conversion already happened globally in P9 — new per-cell loops use
   `include/frehg2/core/ParallelFor.hpp` from the start and must stay bit-identical CPU/OpenMP.
6. Add `tests/<module>/test_*.cpp` with numerical assertions (≥ 2 parameter sets, ≥ 1 edge case,
   ≥ 1 failure-mode test — CQ6). Register them in `tests/CMakeLists.txt`.

---

## 7. Provenance & git

Build provenance (`git_sha`, `git_dirty`) is captured automatically via the CMake-generated
`frehg2/io/build_info.hpp`. Until `git init` is run it reports `unknown`; once the repo is under
git, every output file and checkpoint records the resolved git SHA and config SHA-256 with no
further action.
