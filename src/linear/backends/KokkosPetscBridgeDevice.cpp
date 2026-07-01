// P10 — Option B (GPU-native) Kokkos -> PETSc COO assembly. Tasks 10.3.1 / 10.3.3.
//
// This translation unit is compiled ONLY on a Kokkos GPU backend (CUDA / HIP / SYCL); the
// whole body is wrapped in `#if FREHG2_GPU_ASSEMBLY`, so on the macOS OpenMP/Serial build it
// is an empty object file and the host-staged Option A in KokkosPetscBridge.cpp is used.
//
// INVARIANT (enforced by the `check_gpu_coo_assembly` static check): this file contains no
// host-staging matrix insertion call — the GPU assembly path goes exclusively through the COO
// insertion routine, the only PETSc entry point that accepts device-resident value pointers.
#include "linear/backends/KokkosPetscBridge.hpp"

#if FREHG2_GPU_ASSEMBLY

#include <stdexcept>
#include <vector>

namespace frehg2 {

namespace {
// The device value pointer is handed to PETSc verbatim, so PetscScalar must match `real`.
static_assert(sizeof(PetscScalar) == sizeof(real),
              "Frehg2 requires a real-scalar PETSc build (PetscScalar == double)");

void checkPetscDev(PetscErrorCode ierr, const char* what) {
  if (ierr != 0) {
    throw std::runtime_error(std::string("PETSc (GPU COO) call failed: ") + what);
  }
}
}  // namespace

void KokkosPetscBridge::setPreallocationCOODevice(Mat A, const IntArr1D& rows,
                                                  const IntArr1D& cols) {
  const size_t nnz = rows.extent(0);
  if (cols.extent(0) != nnz) {
    throw std::runtime_error("KokkosPetscBridge::setPreallocationCOODevice: size mismatch");
  }
  // MatSetPreallocationCOO consumes the (i,j) index arrays on the host; only the per-step
  // values (assembleCOODevice) stay device-resident. Stage the indices once here.
  auto rows_h = Kokkos::create_mirror_view(rows);
  auto cols_h = Kokkos::create_mirror_view(cols);
  Kokkos::deep_copy(rows_h, rows);
  Kokkos::deep_copy(cols_h, cols);

  std::vector<PetscInt> coo_i(nnz);
  std::vector<PetscInt> coo_j(nnz);
  for (size_t e = 0; e < nnz; ++e) {
    coo_i[e] = static_cast<PetscInt>(rows_h(e));
    coo_j[e] = static_cast<PetscInt>(cols_h(e));
  }
  checkPetscDev(
      MatSetPreallocationCOO(A, static_cast<PetscCount>(nnz), coo_i.data(), coo_j.data()),
      "MatSetPreallocationCOO");
}

void KokkosPetscBridge::assembleCOODevice(Mat A, const RealArr1D& vals) {
  // vals lives in device memory; reinterpret its device pointer as PetscScalar* and let
  // PETSc (Kokkos-aware MATAIJKOKKOS) consume it directly. Duplicate (i,j) entries from the
  // preallocated pattern are summed by PETSc. No host staging, no MatSetValues.
  const PetscScalar* v = reinterpret_cast<const PetscScalar*>(vals.data());
  checkPetscDev(MatSetValuesCOO(A, v, INSERT_VALUES), "MatSetValuesCOO");
}

}  // namespace frehg2

#endif  // FREHG2_GPU_ASSEMBLY
