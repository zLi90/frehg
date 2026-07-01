// Kokkos -> PETSc assembly bridge (P2.6). Lives INSIDE the PETSc backend; physics code
// never sees it (it only calls SparseSystem::addRow/addCOO).
//
// Two paths: Option A stages coefficients through host memory and calls MatSetValues
// (works for every PETSc Mat type, CPU and GPU). Option B (MatSetValuesCOO, device
// pointers) is GPU-native and requires PETSc >= 3.17 built with Kokkos; it is wired in
// P10. CRITICAL: MatSetValues must NOT be given device pointers on CPU matrix types, so
// the bridge always stages through host memory here.
#ifndef FREHG2_LINEAR_BACKENDS_KOKKOS_PETSC_BRIDGE_HPP
#define FREHG2_LINEAR_BACKENDS_KOKKOS_PETSC_BRIDGE_HPP

#include <petscmat.h>

#include "frehg2/core/define.hpp"
#include "linear/backends/GpuAssembly.hpp"

namespace frehg2 {

class KokkosPetscBridge {
 public:
  // Option A: insert one row from host-resident global indices/values via MatSetValues.
  static void insertRowHost(Mat A, int global_row, int n_cols, const int* global_cols_host,
                            const real* vals_host, InsertMode mode);

  // Host-staged COO assembly (CPU path / Option A for whole-matrix COO). Each (row,col)
  // entry is inserted via MatSetValues. Views must be HostSpace.
  static void assembleCOOHost(Mat A, const IntArr1DHost& rows, const IntArr1DHost& cols,
                              const RealArr1DHost& vals, InsertMode mode);

  // Convenience: host pointer of a HostSpace view (used by tests to feed MatSetValues).
  static real* getHostPtr(RealArr1DHost& v) { return v.data(); }

#if FREHG2_GPU_ASSEMBLY
  // ---- Option B (GPU-native, MatSetValuesCOO) — Tasks 10.3.1 / 10.3.3 ----
  // These are compiled ONLY on a Kokkos GPU backend (see GpuAssembly.hpp). They consume
  // device-resident Kokkos Views directly; PETSc reads the device pointers (no host stage).
  // The implementation lives in KokkosPetscBridgeDevice.cpp, which contains NO MatSetValues
  // call (verified by the check_gpu_coo_assembly static check) — only MatSetValuesCOO.

  // Fix the COO sparsity pattern once from device-resident global (row,col) index views.
  // Must be called before the first assembleCOODevice(); PETSc sums duplicate (i,j).
  static void setPreallocationCOODevice(Mat A, const IntArr1D& rows, const IntArr1D& cols);

  // Insert device-resident COO values (same length/order as the pattern) via MatSetValuesCOO.
  static void assembleCOODevice(Mat A, const RealArr1D& vals);
#endif  // FREHG2_GPU_ASSEMBLY
};

}  // namespace frehg2

#endif  // FREHG2_LINEAR_BACKENDS_KOKKOS_PETSC_BRIDGE_HPP
