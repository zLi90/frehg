// P2.5 acceptance: the assembly kernel emits COO triplets through the SparseSystem
// interface with NO PETSc/solver-library types in its signature or body. The function
// below (assembleSpdSystem) is exactly the shape a physics module uses; it includes only
// the backend-agnostic interface. The concrete PetscLinearSolver is instantiated only in
// the test driver to solve the assembled system. When a second backend is added (P2B),
// the SAME assembleSpdSystem feeds it unchanged and must match within L2 < 1e-10.
#include <cmath>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <petscsys.h>

#include "frehg2_test.hpp"
#include "frehg2/linear/LinearSolver.hpp"
#include "frehg2/linear/SparseSystem.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"

using namespace frehg2;

namespace {
real xExact(int grow) { return 0.5 + 0.001 * static_cast<real>(grow); }

// *** Backend-agnostic assembly: only DecompBase + SparseSystem appear. No PETSc. ***
void assembleSpdSystem(const DecompBase& dd, SparseSystem& sys, const RealArr1D& b_owned) {
  sys.zeroEntries();
  sys.beginAssembly();
  const int n_local = dd.ownedRowCount();
  auto bh = Kokkos::create_mirror_view(b_owned);
  int cols[kMaxStencil];
  for (int L = 0; L < n_local; ++L) {
    int grow = 0, ncols = 0;
    dd.ownedRowStencil(L, grow, cols, ncols);
    real vals[kMaxStencil];
    vals[0] = 4.0;
    real b = 4.0 * xExact(grow);
    for (int c = 1; c < ncols; ++c) {
      vals[c] = -1.0;
      b += -1.0 * xExact(cols[c]);
    }
    sys.addRow(grow, ncols, cols, vals);
    bh(static_cast<size_t>(L)) = b;
  }
  Kokkos::deep_copy(b_owned, bh);
  sys.endAssembly();
}
}  // namespace

TEST_CASE("backend-agnostic assembly solved through PETSc backend") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int mpi_nx = (size >= 2) ? 2 : 1;
  int mpi_ny = 1;
  if (mpi_nx * mpi_ny != size) { mpi_nx = size; }
  MpiComm mc(9, 7, mpi_nx, mpi_ny);
  Decomp2D dd(mc);

  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = (size == 1) ? "jacobi" : "gamg";
  cfg.rtol = 1e-13;
  PetscLinearSolver solver(cfg);
  auto sys = solver.createSystem(dd);

  RealArr1D b("b", dd.ownedRowCount());
  RealArr1D x("x", dd.ownedRowCount());
  assembleSpdSystem(dd, *sys, b);  // <-- no PETSc types crossed this boundary
  solver.solve(*sys, b, x);

  auto xh = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(xh, x);
  const int rstart = dd.ownershipRange().first;
  real local = 0.0;
  for (int L = 0; L < dd.ownedRowCount(); ++L) {
    const real e = xh(static_cast<size_t>(L)) - xExact(rstart + L);
    local += e * e;
  }
  real global = local;
  MPI_Allreduce(MPI_IN_PLACE, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  REQUIRE(std::sqrt(global) < 1e-8);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rc = frehg2test::runAll();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  MPI_Finalize();
  Kokkos::finalize();
  return grc;
}
