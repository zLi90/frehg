// Default (and only currently sanctioned) production linear-solver backend (P2.5.4).
//
// Builds a MATMPIAIJ (or MATAIJKOKKOS when FREHG2_USE_GPU_MAT) matrix sized and
// preallocated from the model-owned DomainDecomposition (no DMDA, no sequential-only
// matrix creation, no self-communicator). Serial == one-rank MPIAIJ. PETSc types appear
// here because this file lives under src/linear/backends/ (the only place the seam check
// permits them).
#ifndef FREHG2_LINEAR_BACKENDS_PETSC_LINEAR_SOLVER_HPP
#define FREHG2_LINEAR_BACKENDS_PETSC_LINEAR_SOLVER_HPP

#include <petscksp.h>

#include <memory>

#include "frehg2/core/define.hpp"
#include "frehg2/linear/LinearSolver.hpp"
#include "frehg2/linear/SolverConfig.hpp"
#include "frehg2/linear/SparseSystem.hpp"
#include "linear/backends/GpuAssembly.hpp"

namespace frehg2 {

class DecompBase;

// PETSc-backed sparse system (MATMPIAIJ + two MPI Vecs for rhs/solution).
class PetscSparseSystem : public SparseSystem {
 public:
  explicit PetscSparseSystem(const DecompBase& dd);
  ~PetscSparseSystem() override;

  PetscSparseSystem(const PetscSparseSystem&) = delete;
  PetscSparseSystem& operator=(const PetscSparseSystem&) = delete;

  void preallocate(const DecompBase& dd) override;
  void beginAssembly() override;
  void endAssembly() override;
  void addRow(int global_row, int ncols, const int* global_cols, const real* vals) override;
  void addCOO(const IntArr1D& rows, const IntArr1D& cols, const RealArr1D& vals) override;
  void setRhs(const RealArr1D& b_owned) override;
  void getSolution(RealArr1D& x_owned) override;
  void zeroEntries() override;

  // Backend-internal accessors (used by PetscLinearSolver::solve).
  Mat mat() const { return A_; }
  Vec rhsVec() const { return b_; }
  Vec solVec() const { return x_; }
  PetscInt localRows() const { return n_local_; }
  PetscInt rowStart() const { return row_start_; }

 private:
  Mat A_ = nullptr;
  Vec b_ = nullptr;
  Vec x_ = nullptr;
  PetscInt n_local_ = 0;
  PetscInt n_global_ = 0;
  PetscInt row_start_ = 0;
#if FREHG2_GPU_ASSEMBLY
  // Option B keeps the COO sparsity pattern resident so MatSetPreallocationCOO runs once and
  // every step only ships device values via MatSetValuesCOO (set on the first addCOO call).
  bool coo_preallocated_ = false;
#endif
};

class PetscLinearSolver : public LinearSolver {
 public:
  explicit PetscLinearSolver(SolverConfig cfg = SolverConfig());
  ~PetscLinearSolver() override;

  PetscLinearSolver(const PetscLinearSolver&) = delete;
  PetscLinearSolver& operator=(const PetscLinearSolver&) = delete;

  std::unique_ptr<SparseSystem> createSystem(const DecompBase& dd) override;
  void setup(SparseSystem& A) override;
  void solve(SparseSystem& A, const RealArr1D& b, RealArr1D& x) override;
  int getIterationCount() const override { return last_iters_; }
  real getResidualNorm() const override { return last_resnorm_; }

  const SolverConfig& config() const { return cfg_; }

 private:
  void ensureKsp();

  SolverConfig cfg_;
  KSP ksp_ = nullptr;
  int last_iters_ = 0;
  real last_resnorm_ = 0.0;
};

}  // namespace frehg2

#endif  // FREHG2_LINEAR_BACKENDS_PETSC_LINEAR_SOLVER_HPP
