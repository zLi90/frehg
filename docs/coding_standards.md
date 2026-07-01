# Frehg2 Coding Standards

These rules are derived from [`.cursorrules`](../.cursorrules) (the authoritative engineering
policy) and the realized code base. A change that violates any of them fails review. The CMake
seam-check targets and CI enforce the structural rules automatically.

---

## 1. Language & formatting

- **C++20** (`CMAKE_CXX_STANDARD 20`, no extensions). Do not use C++23-only features without
  explicit approval. Kokkos is built with `Kokkos_CXX_STANDARD=20`.
- **Formatting is mechanical** — `.clang-format` at the repo root is authoritative (LLVM base,
  2-space indent, 100-column, attached braces, non-indented namespaces, left-aligned pointers):

  ```bash
  clang-format -i $(git ls-files '*.hpp' '*.cpp')   # format
  clang-format --dry-run --Werror $(git ls-files '*.hpp' '*.cpp')   # CI check
  ```

- **Warnings are errors.** `FREHG2_STRICT_WARNINGS=ON` compiles with
  `-Wall -Wextra -Werror -Wno-unused-parameter`. If a warning is truly unavoidable, suppress it
  narrowly with `#pragma GCC diagnostic` and a comment — do not disable `-Werror`.

---

## 2. Naming

| Kind | Convention | Example |
|------|------------|---------|
| Function / method | `lowerCamelCase` | `computeExchangeRate` |
| Class / struct | `UpperCamelCase` | `SweSolver` |
| Member variable | trailing underscore | `state_` |
| Constant | `UPPER_SNAKE_CASE` | `MIN_DEPTH` |
| Floating-point type | the `real` typedef (default `double`) | `real eta;` |

**Variable convention is locked (DP4):** surface `eta, u, v`; subsurface `pressure_head`; solute
`conc`. The surface matrix is `A · eta_new = b`. The SERGHEI names `h, hu, hv` (and the ambiguous
`H`, `q_x/q_y/q_z`, `K_harmonic`) are **forbidden** in new code.

---

## 3. File organization (CQ3)

- Public API in `include/frehg2/<module>/*.hpp`; implementation in `src/<module>/*.cpp`.
- **Every `.cpp` has a corresponding `.hpp`** (except `main.cpp`). No header-only classes — the
  file layout is enforced by a CMake glob; a missing `.cpp` is a build error.
- `#include <Kokkos_Core.hpp>` at the top of any file using Kokkos. `#include <petsc.h>` (and any
  solver-library header) is allowed **only** inside `src/linear/backends/`.

---

## 4. Code-quality rules (CQ1–CQ6)

- **CQ1 — No dead code.** No commented-out code, no `#if 0`, no unused variables. Delete it; do not
  keep "might need this" branches.
- **CQ2 — Header guards, not `#pragma once`.** Use `#ifndef FREHG2_<MODULE>_<NAME>_HPP` guards.
- **CQ3 — Separate `.hpp` / `.cpp`** (see §3).
- **CQ4 — No globals, no singletons.** All state lives in `State`, `GwState`, or objects owned by
  `Orchestrator`. No `static T instance;`. Test isolation depends on this.
- **CQ5 — No stubs in the production path.** A production function is implemented or absent. No
  `return 0;` / `return {};` / `throw std::runtime_error("TODO")` and no `// TODO`/`// FIXME` in
  `src/`. Stubs are allowed only in `tests/` for mocks.
- **CQ6 — Tests cover failure modes.** Every error path has an explicit test asserting the error
  (CFL refusal, malformed YAML, missing IC file, restart mismatch, schema-version refusal, …).

---

## 5. Memory safety & error handling

- **Arrays:** `Kokkos::View` only — no `new[]`/`delete[]`/`malloc`/`free`, no raw owning pointers,
  no `std::vector<real>` for solver state (not GPU-capable / no dual space).
- **PETSc objects:** RAII wrappers / custom deleters (in the backend). HDF5 ids: `h5::Guard`.
- **Errors:** `PetscErrorCode` + `CHKERRQ` for PETSc; C++ exceptions on the host (never inside a
  Kokkos lambda); `Kokkos::abort()` for device-side errors. Messages say *what* failed, *where*,
  and *why*.

---

## 6. Kokkos / GPU rules (DP3)

- Per-cell loops go through `include/frehg2/core/ParallelFor.hpp` wrappers (`parallelForRange`/
  `Surface`/`Volume`). Reductions use `Kokkos::Max`/`Min` (order-independent ⇒ bit-identical).
- Inside `KOKKOS_LAMBDA` / device-callable helpers: **no** `std::cout`/`printf`/`std::vector`/file
  I/O, **no** `MatSetValues` (PETSc is not GPU-callable), `std::` math → `Kokkos::` math. Device
  helpers are `KOKKOS_INLINE_FUNCTION`.
- Physics kernels emit COO triplets to `SparseSystem::addRow`(host)/`addCOO`(device); the
  Kokkos→PETSc staging lives **only** in `src/linear/backends/`.
- Loops kept deliberately sequential (Gauss–Seidel sweeps, scatter limiters, prefix sums, COO
  assembly) are catalogued in [`local_loops_audit.md`](local_loops_audit.md) with the reason.

---

## 7. Linear-algebra & parallel rules (DP2)

- All global solves go through the `LinearSolver`/`SparseSystem` interface. Physics never sees
  `Mat`/`Vec`/`KSP`/`DM`/`Tpetra::*`.
- **No PETSc `DMDA`/`DMCreateMatrix`/`DMDACreate*`.** Matrices are `MATMPIAIJ` (or `MATAIJKOKKOS`)
  preallocated from `DomainDecomposition`. **No `MatCreateSeqAIJ`/`PETSC_COMM_SELF`** — serial is
  one-rank MPIAIJ. **No hand-rolled Krylov solvers** — PETSc KSP is the only sanctioned backend.
- A kernel reading neighbor cells must call the model-owned halo exchange first
  (`MpiComm::haloExchange*`); reading a ghost without exchanging is a data race.
- MPI rank-count equivalence (1/2/4 ranks, same executable) must agree to `L2 < 1e-10`.

These are enforced at build time by `check_no_seqaij`, `check_solver_seam`, and
`check_gpu_coo_assembly` (run as dependencies of the `frehg2` target).

---

## 8. Testing standards

Every test asserts a numerical bound (`REQUIRE(error < tol)`); "it compiles" / "no crash" is not a
pass. Prefer ≥ 2 parameter sets and ≥ 1 edge case (dry cells, zero flux, NODATA). Default
tolerances: matrix/RHS regression `< 1e-9`; one-step solver regression `< 1e-8`; full benchmark
time series `< 1e-6`; MPI rank equivalence `< 1e-10`. Do not loosen a benchmark gate or defer
required behavior to a later phase to declare completion.
