#include "linear/backends/PetscLinearSolver.hpp"

#include <array>
#include <stdexcept>
#include <vector>

#include "linear/DomainDecomposition.hpp"
#include "linear/backends/KokkosPetscBridge.hpp"

namespace frehg2 {

namespace {
#ifdef FREHG2_USE_GPU_MAT
constexpr const char* kMatType = MATAIJKOKKOS;
constexpr const char* kVecType = VECKOKKOS;
#else
constexpr const char* kMatType = MATMPIAIJ;  // serial == one-rank MPIAIJ (never seqaij)
constexpr const char* kVecType = VECMPI;
#endif

void checkPetsc(PetscErrorCode ierr, const char* what) {
  if (ierr != 0) {
    throw std::runtime_error(std::string("PETSc call failed: ") + what);
  }
}
}  // namespace

// ================================ PetscSparseSystem ================================

PetscSparseSystem::PetscSparseSystem(const DecompBase& dd) {
  n_local_ = dd.ownedRowCount();
  n_global_ = dd.globalRowCount();
  row_start_ = dd.ownershipRange().first;

  checkPetsc(MatCreate(PETSC_COMM_WORLD, &A_), "MatCreate");
  checkPetsc(MatSetSizes(A_, n_local_, n_local_, n_global_, n_global_), "MatSetSizes");
  checkPetsc(MatSetType(A_, kMatType), "MatSetType");

  checkPetsc(VecCreate(PETSC_COMM_WORLD, &b_), "VecCreate");
  checkPetsc(VecSetSizes(b_, n_local_, n_global_), "VecSetSizes");
  checkPetsc(VecSetType(b_, kVecType), "VecSetType");
  checkPetsc(VecDuplicate(b_, &x_), "VecDuplicate");
}

PetscSparseSystem::~PetscSparseSystem() {
  MatDestroy(&A_);
  VecDestroy(&b_);
  VecDestroy(&x_);
}

void PetscSparseSystem::preallocate(const DecompBase& dd) {
  std::vector<PetscInt> d_nnz(static_cast<size_t>(n_local_), 0);
  std::vector<PetscInt> o_nnz(static_cast<size_t>(n_local_), 0);
  const PetscInt rlo = row_start_;
  const PetscInt rhi = row_start_ + n_local_;

  std::array<int, kMaxStencil> cols{};
  for (int L = 0; L < n_local_; ++L) {
    int grow = 0;
    int ncols = 0;
    dd.ownedRowStencil(L, grow, cols.data(), ncols);
    for (int c = 0; c < ncols; ++c) {
      const PetscInt gc = cols[static_cast<size_t>(c)];
      if (gc >= rlo && gc < rhi) {
        d_nnz[static_cast<size_t>(L)]++;
      } else {
        o_nnz[static_cast<size_t>(L)]++;
      }
    }
  }

  // Works for MATMPIAIJ and the AIJ-derived MATAIJKOKKOS, on one or many ranks.
  checkPetsc(MatMPIAIJSetPreallocation(A_, 0, d_nnz.data(), 0, o_nnz.data()),
             "MatMPIAIJSetPreallocation");
  // One-rank MPIAIJ keeps everything in the diagonal block; also set the seq path so a
  // single-rank communicator is fully preallocated.
  checkPetsc(MatSeqAIJSetPreallocation(A_, 0, d_nnz.data()), "MatSeqAIJSetPreallocation");
  checkPetsc(MatSetUp(A_), "MatSetUp");
}

void PetscSparseSystem::beginAssembly() {
  // No state needed; addRow/addCOO insert directly. Kept for interface symmetry.
}

void PetscSparseSystem::endAssembly() {
  checkPetsc(MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY), "MatAssemblyBegin");
  checkPetsc(MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY), "MatAssemblyEnd");
}

void PetscSparseSystem::addRow(int global_row, int ncols, const int* global_cols,
                               const real* vals) {
  KokkosPetscBridge::insertRowHost(A_, global_row, ncols, global_cols, vals, ADD_VALUES);
}

void PetscSparseSystem::addCOO(const IntArr1D& rows, const IntArr1D& cols,
                               const RealArr1D& vals) {
#if FREHG2_GPU_ASSEMBLY
  // Option B (GPU-native): fix the (i,j) pattern once, then ship device-resident values via
  // MatSetValuesCOO every step — zero host staging. The physics emits the same pattern each
  // step, so the preallocation is cached after the first call.
  if (!coo_preallocated_) {
    KokkosPetscBridge::setPreallocationCOODevice(A_, rows, cols);
    coo_preallocated_ = true;
  }
  KokkosPetscBridge::assembleCOODevice(A_, vals);
#else
  // Option A (default, CPU/OpenMP): stage device views to host and insert per entry. This is
  // bit-identical to the addRow path used by every CPU gate.
  auto rh = Kokkos::create_mirror_view(rows);
  auto ch = Kokkos::create_mirror_view(cols);
  auto vh = Kokkos::create_mirror_view(vals);
  Kokkos::deep_copy(rh, rows);
  Kokkos::deep_copy(ch, cols);
  Kokkos::deep_copy(vh, vals);
  KokkosPetscBridge::assembleCOOHost(A_, rh, ch, vh, ADD_VALUES);
#endif
}

void PetscSparseSystem::setRhs(const RealArr1D& b_owned) {
#if FREHG2_GPU_ASSEMBLY
  // Device Vec bridge (Task 10.3.4): copy device->device into the VECKOKKOS array, keeping
  // the RHS resident on the GPU (no host round-trip). PETSc returns the device pointer for a
  // Kokkos vector via VecGetArrayAndMemType.
  PetscScalar* arr = nullptr;
  PetscMemType mtype;
  checkPetsc(VecGetArrayAndMemType(b_, &arr, &mtype), "VecGetArrayAndMemType(rhs)");
  Kokkos::View<real*, Kokkos::LayoutRight, DeviceSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      dst(reinterpret_cast<real*>(arr), static_cast<size_t>(n_local_));
  Kokkos::deep_copy(dst, Kokkos::subview(b_owned, Kokkos::make_pair(0, static_cast<int>(n_local_))));
  checkPetsc(VecRestoreArrayAndMemType(b_, &arr), "VecRestoreArrayAndMemType(rhs)");
#else
  RealArr1DHost b_host = b_owned;
  if (b_host.data() != b_owned.data()) Kokkos::deep_copy(b_host, b_owned);
  PetscScalar* arr = nullptr;
  checkPetsc(VecGetArray(b_, &arr), "VecGetArray(rhs)");
  for (PetscInt L = 0; L < n_local_; ++L) {
    arr[L] = static_cast<PetscScalar>(b_host(static_cast<size_t>(L)));
  }
  checkPetsc(VecRestoreArray(b_, &arr), "VecRestoreArray(rhs)");
#endif
}

void PetscSparseSystem::getSolution(RealArr1D& x_owned) {
#if FREHG2_GPU_ASSEMBLY
  const PetscScalar* arr = nullptr;
  PetscMemType mtype;
  checkPetsc(VecGetArrayReadAndMemType(x_, &arr, &mtype), "VecGetArrayReadAndMemType(sol)");
  Kokkos::View<const real*, Kokkos::LayoutRight, DeviceSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      src(reinterpret_cast<const real*>(arr), static_cast<size_t>(n_local_));
  Kokkos::deep_copy(Kokkos::subview(x_owned, Kokkos::make_pair(0, static_cast<int>(n_local_))),
                    src);
  checkPetsc(VecRestoreArrayReadAndMemType(x_, &arr), "VecRestoreArrayReadAndMemType(sol)");
#else
  RealArr1DHost x_host = x_owned;
  const PetscScalar* arr = nullptr;
  checkPetsc(VecGetArrayRead(x_, &arr), "VecGetArrayRead(sol)");
  for (PetscInt L = 0; L < n_local_; ++L) {
    x_host(static_cast<size_t>(L)) = static_cast<real>(PetscRealPart(arr[L]));
  }
  checkPetsc(VecRestoreArrayRead(x_, &arr), "VecRestoreArrayRead(sol)");
  if (x_host.data() != x_owned.data()) Kokkos::deep_copy(x_owned, x_host);
#endif
}

void PetscSparseSystem::zeroEntries() { checkPetsc(MatZeroEntries(A_), "MatZeroEntries"); }

// ================================ PetscLinearSolver ================================

PetscLinearSolver::PetscLinearSolver(SolverConfig cfg) : cfg_(std::move(cfg)) {}

PetscLinearSolver::~PetscLinearSolver() { KSPDestroy(&ksp_); }

void PetscLinearSolver::ensureKsp() {
  if (ksp_ == nullptr) {
    checkPetsc(KSPCreate(PETSC_COMM_WORLD, &ksp_), "KSPCreate");
  }
}

std::unique_ptr<SparseSystem> PetscLinearSolver::createSystem(const DecompBase& dd) {
  auto sys = std::make_unique<PetscSparseSystem>(dd);
  sys->preallocate(dd);
  return sys;
}

void PetscLinearSolver::setup(SparseSystem& A) {
  ensureKsp();
  auto& sys = dynamic_cast<PetscSparseSystem&>(A);
  checkPetsc(KSPSetOperators(ksp_, sys.mat(), sys.mat()), "KSPSetOperators");
  checkPetsc(KSPSetType(ksp_, cfg_.ksp_type.c_str()), "KSPSetType");
  PC pc = nullptr;
  checkPetsc(KSPGetPC(ksp_, &pc), "KSPGetPC");
  checkPetsc(PCSetType(pc, cfg_.pc_type.c_str()), "PCSetType");
  checkPetsc(KSPSetTolerances(ksp_, cfg_.rtol, cfg_.atol, PETSC_DEFAULT, cfg_.max_it),
             "KSPSetTolerances");
  // Allow command-line/options-db overrides (e.g. -ksp_type, -pc_type) on top of cfg.
  checkPetsc(KSPSetFromOptions(ksp_), "KSPSetFromOptions");
}

void PetscLinearSolver::solve(SparseSystem& A, const RealArr1D& b, RealArr1D& x) {
  auto& sys = dynamic_cast<PetscSparseSystem&>(A);
  setup(A);  // idempotent; binds current operator and applies config
  sys.setRhs(b);
  checkPetsc(KSPSolve(ksp_, sys.rhsVec(), sys.solVec()), "KSPSolve");

  KSPConvergedReason reason;
  checkPetsc(KSPGetConvergedReason(ksp_, &reason), "KSPGetConvergedReason");
  if (reason < 0) {
    throw std::runtime_error("PetscLinearSolver::solve: KSP diverged (reason " +
                             std::to_string(static_cast<int>(reason)) + ")");
  }
  PetscInt its = 0;
  PetscReal rnorm = 0.0;
  checkPetsc(KSPGetIterationNumber(ksp_, &its), "KSPGetIterationNumber");
  checkPetsc(KSPGetResidualNorm(ksp_, &rnorm), "KSPGetResidualNorm");
  last_iters_ = static_cast<int>(its);
  last_resnorm_ = static_cast<real>(rnorm);

  sys.getSolution(x);
}

}  // namespace frehg2
