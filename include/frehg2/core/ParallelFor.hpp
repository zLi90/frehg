// Frehg2 on-node parallel-loop helpers (P9).
//
// Phase 9 is the single, last pass that converts the per-cell *local* kernels (SWE/RE flux &
// geometry, source/sink, coupling exchange, solute advection/diffusion fills, and I/O packing)
// from plain `for` loops to `Kokkos::parallel_for` / `Kokkos::parallel_reduce`. It introduces
// NO new numerical behaviour: each kernel still performs the identical per-cell arithmetic, so
// on the host (OpenMP/Serial) backend the results are bit-identical to the P4/P5 serial
// references (gate). Only loops whose iterations write *disjoint* cells (or are pure
// max/min/sum reductions) are converted here; order-dependent loops (PETSc COO assembly,
// in-place neighbour writes, MPI halo packing, 1-D boundary-edge ghost fills) stay sequential.
//
// The execution space is the default HOST space because Frehg2 state currently lives in
// `RealArr1DHost` (HostSpace) views. The P10 GPU bridge migrates the state to device views and
// flips `LoopExec` to `Kokkos::DefaultExecutionSpace`; the kernel bodies below are already
// written as `KOKKOS_LAMBDA`s capturing only value types (Grid, scalars, View handles), so that
// migration needs no kernel rewrite.
#ifndef FREHG2_CORE_PARALLEL_FOR_HPP
#define FREHG2_CORE_PARALLEL_FOR_HPP

#include <Kokkos_Core.hpp>

namespace frehg2 {

// Host execution space for the local per-cell kernels (OpenMP on macOS; Serial fallback).
using LoopExec = Kokkos::DefaultHostExecutionSpace;

// Flat range [0, n): used for whole-array element-wise updates (state copies, resets).
template <class Functor>
inline void parallelForRange(const char* name, int n, const Functor& f) {
  Kokkos::parallel_for(name, Kokkos::RangePolicy<LoopExec>(0, n), f);
}

// 2-D surface tile [0,nx) x [0,ny); the lambda is called as f(i, j).
template <class Functor>
inline void parallelForSurface(const char* name, int nx, int ny, const Functor& f) {
  using Policy = Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<2>>;
  Kokkos::parallel_for(name, Policy({0, 0}, {nx, ny}), f);
}

// 3-D subsurface tile [0,nx) x [0,ny) x [0,nz); the lambda is called as f(i, j, k).
template <class Functor>
inline void parallelForVolume(const char* name, int nx, int ny, int nz, const Functor& f) {
  using Policy = Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<3>>;
  Kokkos::parallel_for(name, Policy({0, 0, 0}, {nx, ny, nz}), f);
}

}  // namespace frehg2

#endif  // FREHG2_CORE_PARALLEL_FOR_HPP
