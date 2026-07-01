// Backend-agnostic linear-solver interface (P2.5.2).
//
// Physics code holds a LinearSolver& and a SparseSystem; the concrete backend (PETSc by
// default; optionally Trilinos/Ginkgo/KokkosKernels per Phase 2B) is chosen at
// construction. No solver-library type leaks through this interface.
#ifndef FREHG2_LINEAR_LINEAR_SOLVER_HPP
#define FREHG2_LINEAR_LINEAR_SOLVER_HPP

#include <memory>

#include "frehg2/core/define.hpp"
#include "frehg2/linear/SparseSystem.hpp"

namespace frehg2 {

class DecompBase;

class LinearSolver {
 public:
  virtual ~LinearSolver() = default;

  // Create a backend-native sparse system sized/structured for the decomposition.
  virtual std::unique_ptr<SparseSystem> createSystem(const DecompBase& dd) = 0;

  // Bind the (assembled) operator and configure the Krylov method / preconditioner.
  virtual void setup(SparseSystem& A) = 0;

  // Solve A x = b for the owned unknowns (b/x in owned contiguous ordering).
  virtual void solve(SparseSystem& A, const RealArr1D& b, RealArr1D& x) = 0;

  virtual int getIterationCount() const = 0;
  virtual real getResidualNorm() const = 0;
};

}  // namespace frehg2

#endif  // FREHG2_LINEAR_LINEAR_SOLVER_HPP
