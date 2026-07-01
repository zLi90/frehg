#include "linear/backends/KokkosPetscBridge.hpp"

#include <array>
#include <stdexcept>

#include "linear/DomainDecomposition.hpp"  // for kMaxStencil

namespace frehg2 {

namespace {
// PETSc must be a real-scalar build for the direct real* -> PetscScalar* reinterpretation.
static_assert(sizeof(PetscScalar) == sizeof(real),
              "Frehg2 requires a real-scalar PETSc build (PetscScalar == double)");
}  // namespace

void KokkosPetscBridge::insertRowHost(Mat A, int global_row, int n_cols,
                                      const int* global_cols_host, const real* vals_host,
                                      InsertMode mode) {
  // Convert int -> PetscInt on a small stack buffer (n_cols <= STAR stencil width).
  std::array<PetscInt, kMaxStencil> cols{};
  if (n_cols > static_cast<int>(cols.size())) {
    throw std::runtime_error("KokkosPetscBridge::insertRowHost: too many columns");
  }
  for (int c = 0; c < n_cols; ++c) {
    cols[static_cast<size_t>(c)] = static_cast<PetscInt>(global_cols_host[c]);
  }
  PetscInt row = static_cast<PetscInt>(global_row);
  const PetscScalar* vals = reinterpret_cast<const PetscScalar*>(vals_host);
  PetscErrorCode ierr = MatSetValues(A, 1, &row, n_cols, cols.data(), vals, mode);
  if (ierr != 0) {
    throw std::runtime_error("KokkosPetscBridge::insertRowHost: MatSetValues failed");
  }
}

void KokkosPetscBridge::assembleCOOHost(Mat A, const IntArr1DHost& rows,
                                        const IntArr1DHost& cols, const RealArr1DHost& vals,
                                        InsertMode mode) {
  const size_t nnz = rows.extent(0);
  if (cols.extent(0) != nnz || vals.extent(0) != nnz) {
    throw std::runtime_error("KokkosPetscBridge::assembleCOOHost: size mismatch");
  }
  for (size_t e = 0; e < nnz; ++e) {
    PetscInt r = static_cast<PetscInt>(rows(e));
    PetscInt c = static_cast<PetscInt>(cols(e));
    PetscScalar v = static_cast<PetscScalar>(vals(e));
    PetscErrorCode ierr = MatSetValues(A, 1, &r, 1, &c, &v, mode);
    if (ierr != 0) {
      throw std::runtime_error("KokkosPetscBridge::assembleCOOHost: MatSetValues failed");
    }
  }
}

}  // namespace frehg2
