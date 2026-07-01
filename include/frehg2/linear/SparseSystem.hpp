// Backend-agnostic sparse linear system interface (P2.5.1).
//
// This is the ONLY assembly surface the physics kernels (P4/P5/P8) target. They emit COO
// triplets (global_row, global_cols, vals) via addRow/addCOO; the concrete backend turns
// them into its native matrix. No Mat/Vec/KSP/Tpetra type appears here.
#ifndef FREHG2_LINEAR_SPARSE_SYSTEM_HPP
#define FREHG2_LINEAR_SPARSE_SYSTEM_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

class DecompBase;  // model-owned parallel layout (src/linear/DomainDecomposition.hpp)

class SparseSystem {
 public:
  virtual ~SparseSystem() = default;

  // Set the nonzero pattern from dd.ownedRowStencil over owned rows (the backend computes
  // its own preallocation). Must be called once before assembly.
  virtual void preallocate(const DecompBase& dd) = 0;

  virtual void beginAssembly() = 0;
  virtual void endAssembly() = 0;

  // Host-staged COO add for one matrix row (the only call the CPU assembly kernels use).
  virtual void addRow(int global_row, int ncols, const int* global_cols,
                      const real* vals) = 0;

  // Device-resident COO add (GPU path). Default backends host-stage and reuse addRow.
  virtual void addCOO(const IntArr1D& rows, const IntArr1D& cols,
                      const RealArr1D& vals) = 0;

  // RHS / solution in the owned contiguous ordering (matches DecompBase ownership range).
  virtual void setRhs(const RealArr1D& b_owned) = 0;
  virtual void getSolution(RealArr1D& x_owned) = 0;

  // Reuse the pattern across time steps without reallocating.
  virtual void zeroEntries() = 0;
};

}  // namespace frehg2

#endif  // FREHG2_LINEAR_SPARSE_SYSTEM_HPP
