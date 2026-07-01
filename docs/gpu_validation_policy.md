# GPU Validation Policy (P0.0.6)

Authoritative policy for how GPU capability is developed and validated, given that the
local development machine is macOS (Apple Silicon, no CUDA).

## Principles

1. **GPU-capable code is required throughout.** All per-cell compute is written so it can
   run on device: `Kokkos::View` storage, `KOKKOS_LAMBDA`/`parallel_for` kernels, no
   `std::cout`/`printf`/`std::vector`/file-I/O inside device lambdas, and the
   Kokkos→PETSc bridge uses `MatSetValues` (HostSpace) on CPU and `MatSetValuesCOO`
   (device) on GPU. Matrix/vector types are selected via PETSc options
   (`-mat_type`/`-vec_type`) so the same `LinearSolver` abstraction runs CPU and GPU.

2. **Local macOS gates compile and test CPU / OpenMP / MPI paths only.** The Kokkos
   default execution space on this machine is `Serial`/`OpenMP`. PETSc here is a CPU build
   (3.25.1). All blocking benchmark gates (b1-sw, b2-gw, …) are validated on CPU/OpenMP.

3. **GPU tests are written but skipped locally.** GPU-specific test sources may be added
   and **must** be labeled CTest label `gpu`. They are excluded from the default local
   `ctest` run and are not required to pass on macOS. CI/local builds must not fail
   because a `gpu` test cannot execute.

4. **Phases needing GPU readiness do static checks locally, defer real execution.**
   P9 (Kokkos pass), P10 (GPU backend, incl. GPU async coupling path), P11, P19, P21 may
   enforce GPU-ready code structure via compile-time checks, static greps (no forbidden
   symbols in device code), and CPU-equivalence tests. **Real GPU execution is a deferred
   external validation task** for a future Linux/NVIDIA machine.

5. **CPU/OpenMP completeness gates GPU handoff.** The production CPU/OpenMP model must be
   complete and fully validated (through the b-series gates) before any handoff to Linux
   GPU execution testing. GPU execution must never block completion of the CPU/OpenMP
   model on this machine.

## Operational rules

- Default backend: `PetscLinearSolver` with `-mat_type aij -vec_type standard` on CPU;
  the same code path selects `aijkokkos`/`kokkos` types on a Kokkos-aware PETSc GPU build.
- Forbidden in `src/` (enforced, see `.cursorrules`): `MatCreateSeqAIJ`,
  `MatCreateSeqBAIJ`, `PETSC_COMM_SELF`, `DMDA`/`DMCreateMatrix`-based serial-only paths.
  Serial runs are one MPI rank over `PETSC_COMM_WORLD`.
- A GPU-readiness static check (grep for `std::cout`/`printf`/`std::vector`/`fopen`/
  `MatSetValues` inside `KOKKOS_LAMBDA`, and for forbidden PETSc symbols) is part of the
  P9/P10 acceptance and may be wired as a CTest `gpu`-adjacent lint that runs on CPU.

## Deferred (Linux/NVIDIA) checklist

- Build Kokkos + PETSc with CUDA/HIP; run all `gpu`-labeled tests.
- Validate `MatSetValuesCOO` device assembly equals CPU `MatSetValues` assembly.
- Validate the GPU async coupling path (PetscSubcomm + CUDA/HIP streams, P10.3.7).
- Re-run b-series gates on GPU and confirm CPU/GPU equivalence within solver tolerances.
