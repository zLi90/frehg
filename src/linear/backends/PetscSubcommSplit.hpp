// P10 Task 10.3.7 — GPU realization of the P11 async SW<->GW coupling pipeline.
//
// The CPU async pipeline (src/driver/AsyncPipeline) overlaps the surface and groundwater
// solves on two std::threads but SERIALIZES the actual KSPSolve calls behind a process mutex,
// because PETSc KSP on a SHARED communicator is not thread-safe (and on >1 rank two threads
// issuing collectives on the same comm can deadlock). The fix — needed only once the GPU
// bridge exists — is to give each physics domain its OWN sub-communicator so the two solves
// are independent, and to place each domain's device work on its OWN CUDA/HIP stream so the
// GPU overlaps them. The design is recorded in docs/research_notes/async_gpu.md.
//
// This helper lives in the PETSc backend (the only place solver/PETSc types are permitted) and
// is compiled ONLY on a Kokkos GPU backend (see GpuAssembly.hpp). On the macOS OpenMP build it
// is compiled out entirely, so the production CPU/OpenMP async path (P11) is unaffected and
// remains the numerical reference. Execution is deferred to Linux/NVIDIA per
// docs/gpu_validation_policy.md.
#ifndef FREHG2_LINEAR_BACKENDS_PETSC_SUBCOMM_SPLIT_HPP
#define FREHG2_LINEAR_BACKENDS_PETSC_SUBCOMM_SPLIT_HPP

#include "linear/backends/GpuAssembly.hpp"

#if FREHG2_GPU_ASSEMBLY

#include <petscsys.h>

#include <Kokkos_Core.hpp>

namespace frehg2 {

// Two independent solver groups (surface, groundwater) for the GPU async pipeline: a
// PetscSubcomm split of PETSC_COMM_WORLD plus one stream-backed device execution-space
// instance per group. Each domain's LinearSolver backend builds its Mat/Vec/KSP on its child
// communicator (still MATAIJKOKKOS preallocated from the model-owned DomainDecomposition — no
// DMDA), and launches its kernels + KSPSolve on its own stream, so the two domains run truly
// concurrently. The cross-domain coupling flux is taken on PETSC_COMM_WORLD AFTER both groups
// fence (mirroring the CPU join barrier), so the math is identical to the synchronous path.
class AsyncSolverGroups {
 public:
  enum class Group { Surface = 0, Groundwater = 1 };

  // Split PETSC_COMM_WORLD into two contiguous child communicators and create two
  // stream-backed device exec-space instances (one per domain).
  AsyncSolverGroups();
  ~AsyncSolverGroups();

  AsyncSolverGroups(const AsyncSolverGroups&) = delete;
  AsyncSolverGroups& operator=(const AsyncSolverGroups&) = delete;

  // Child communicator this rank belongs to for the given domain's solver objects.
  MPI_Comm comm(Group g) const;

  // Stream-backed device execution space for the given domain's kernels / assembly.
  const DeviceExec& stream(Group g) const;

  // Barrier equivalent to the CPU pipeline's gw_worker.join(): fence BOTH streams so all
  // device work for the current windows has completed before the cross-domain exchange.
  void fenceAll() const;

 private:
  PetscSubcomm subcomm_ = nullptr;
  DeviceExec stream_surface_;
  DeviceExec stream_groundwater_;
};

}  // namespace frehg2

#endif  // FREHG2_GPU_ASSEMBLY
#endif  // FREHG2_LINEAR_BACKENDS_PETSC_SUBCOMM_SPLIT_HPP
