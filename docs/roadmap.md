# Frehg2 Roadmap (post-1.0)

Frehg2 1.0.0 is the terminal deliverable of [`INTEGRATED_PLAN.md`](../INTEGRATED_PLAN.md) (archived
under [`planning/`](planning/)): a production-grade, GPU-capable, C++20/Kokkos/MPI/PETSc
surface-water — groundwater — solute coupled model with the legacy Frehg algorithm reproduced to
gate tolerance. 1.0 is **not** "feature complete forever". This file records the planned direction.

Items are intentionally scoped, not promised on a date. Each links to a design note where one
exists.

## 1.1 — Solver & discretization strength

- **Tuned algebraic multigrid defaults.** AMG is already reachable at runtime through the
  `LinearSolver`/`SparseSystem` seam (`-pc_type gamg`); 1.1 validates and tunes it as a per-operator
  default on a hypre/MUMPS-enabled build, addressing the P21 KSP hot spot.
  See [`research_notes/multigrid.md`](research_notes/multigrid.md).
- **Adaptive mesh refinement (AMR).** Block-structured patch AMR for localized fronts, integrating
  a mature framework (AMReX) while keeping the physics kernels and the solver seam unchanged.
  See [`research_notes/adaptive_mesh.md`](research_notes/adaptive_mesh.md).
- **Additional constitutive variants.** More Van Genuchten–Mualem variants and alternative
  retention/relative-permeability models behind the existing `SoilParams`/`VanGenuchten` interface.
- **Lateral / layered Richards — DONE (Phase 23).** Horizontal Darcy flux in the
  predictor-corrector RE solver (`groundwater.full_3d: true`) + fully heterogeneous **per-cell**
  soil (3D `SoilMap`), matching legacy `frehg`/SERGHEI (no Picard/Newton). Combined with the new
  partial-width fixed-flux top BC (`groundwater.recharge`, legacy `bctype_GW[5]==2`/`qtop`), this
  unblocked the **full 2-D `b3-kirkland` port** — now a genuine 2-D layered infiltration matching
  the digitized Kirkland `h=0`/`h=-400` contours to ~0.25–0.29 m RMS (review-tier), see
  [`benchmarks/b3-kirkland.md`](benchmarks/b3-kirkland.md). The coupled subsurface leg for
  `b5-vcatchment` (lateral RE under surface flow) remains future work. See `INTEGRATED_PLAN.md`
  Phase 23.
- **Full RE iteration.** 1.0 implements the `iter_solve == 0, use_corrector == 1` path (DP8); 1.1
  generalizes to 2+ corrector passes where a reference justifies it.

## 1.2 — GPU & on-node parallelism

- **GPU execution validation.** Run the `gpu`-labeled tests and the unified validation on a
  Linux/NVIDIA host (`-DFREHG2_ENABLE_CUDA=ON`), the deferred half of P10/P19/P21.
  See [`gpu_validation_policy.md`](gpu_validation_policy.md).
- **Concurrent GPU async coupling.** PetscSubcomm + CUDA-stream split so the SW and GW solves run
  concurrently on the device (Task 10.3.7), promoting the CPU async pipeline (P11) to a real GPU
  speedup. See [`research_notes/async_gpu.md`](research_notes/async_gpu.md).
- **Single-rank multi-GPU** (capacity, not speed), only if profiling on the perf host justifies it
  over the multi-rank/one-GPU path. See
  [`research_notes/single_rank_multi_gpu.md`](research_notes/single_rank_multi_gpu.md).
- **MPI strong-scaling gate.** Parallel efficiency ≥ 0.7 at 16 ranks on `b5-vcatchment` (deferred
  P21 Task 21.3.5) and GPU speedup ≥ 5× at 1M cells (Task 21.3.6), on the HPC/CUDA host.

## 2.0 — SERGHEI-SWE coupling; legacy surface removal

- **Couple with SERGHEI-SWE only.** Replace the legacy-derived surface formulation with a
  SERGHEI-SWE surface, coupling it to the Frehg2 Richards/solute subsurface through the existing
  coupling seam — and **remove the legacy surface convention entirely** (the legacy code archived in
  P20 stops being a reference path).
- **Backend alternatives.** Promote one of the P2B evaluation backends (Trilinos / Ginkgo /
  KokkosKernels) from candidate to a supported production `LinearSolver` backend, if the bake-off
  justifies it. Adopting any requires updating `.cursorrules` and the plan first.

## How this roadmap is maintained

Post-release, each item that did not make a given milestone is tracked as a GitHub issue and
scheduled into the next milestone (P22 Task 22.3.6). This file is the human-readable summary; the
issue tracker is the live status.
