// P2.5 acceptance: solve an SPD system through the LinearSolver/SparseSystem interface.
// Runs on 1 and 2 MPI ranks.
#include <petscksp.h>

#include <cmath>
#include <vector>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"

using namespace frehg2;

namespace {
void procGrid(int size, int& mpi_nx, int& mpi_ny) {
  if (size == 1) { mpi_nx = 1; mpi_ny = 1; }
  else if (size == 2) { mpi_nx = 2; mpi_ny = 1; }
  else if (size == 4) { mpi_nx = 2; mpi_ny = 2; }
  else { mpi_nx = size; mpi_ny = 1; }
}

// Manufactured exact solution as a function of the global row index.
real xExact(int grow) { return 1.5 + std::sin(0.1 * static_cast<real>(grow)); }

// Assemble the 5-point Dirichlet Laplacian (diag 4, in-domain neighbors -1) and the
// consistent RHS b = A x_exact into b_owned.
void assembleLaplacian(Decomp2D& dd, SparseSystem& sys, const RealArr1D& b_owned) {
  sys.zeroEntries();
  sys.beginAssembly();
  const int n_local = dd.ownedRowCount();
  auto bh = Kokkos::create_mirror_view(b_owned);
  int cols[kMaxStencil];
  int ncols = 0;
  for (int L = 0; L < n_local; ++L) {
    int grow = 0;
    dd.ownedRowStencil(L, grow, cols, ncols);
    real vals[kMaxStencil];
    vals[0] = 4.0;  // center
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

real l2ErrorVsExact(Decomp2D& dd, const RealArr1D& x) {
  auto xh = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(xh, x);
  const int rstart = dd.ownershipRange().first;
  const int n_local = dd.ownedRowCount();
  real local = 0.0;
  for (int L = 0; L < n_local; ++L) {
    const real e = xh(static_cast<size_t>(L)) - xExact(rstart + L);
    local += e * e;
  }
  real global = local;
  MPI_Allreduce(MPI_IN_PLACE, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return std::sqrt(global);
}
}  // namespace

TEST_CASE("matrix is MATMPIAIJ and a tight CG solve recovers the exact field (1 rank)") {
  // NOTE: the local PETSc build provides no LU factor for mpiaij (no MUMPS/SuperLU; seq
  // matrices are forbidden by policy). The deterministic serial bring-up path therefore
  // uses CG + Jacobi at a tight tolerance, which recovers the exact solution.
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;
  MpiComm mc(6, 5, 1, 1);
  Decomp2D dd(mc);

  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1e-13;
  PetscLinearSolver solver(cfg);
  auto sys = solver.createSystem(dd);

  // Assert the backend really built an mpiaij matrix.
  auto* psys = dynamic_cast<PetscSparseSystem*>(sys.get());
  REQUIRE(psys != nullptr);
  MatType mtype;
  MatGetType(psys->mat(), &mtype);
  PetscBool is_mpiaij = PETSC_FALSE;
  PetscStrcmp(mtype, MATMPIAIJ, &is_mpiaij);
  REQUIRE(is_mpiaij == PETSC_TRUE);

  RealArr1D b("b", dd.ownedRowCount());
  RealArr1D x("x", dd.ownedRowCount());
  assembleLaplacian(dd, *sys, b);
  solver.setup(*sys);
  solver.solve(*sys, b, x);
  REQUIRE(l2ErrorVsExact(dd, x) < 1e-9);
}

TEST_CASE("cg+gamg: iteration count and residual accessible; solution accurate") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(10, 8, pnx, pny);
  Decomp2D dd(mc);

  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "gamg";
  cfg.rtol = 1e-12;
  PetscLinearSolver solver(cfg);
  auto sys = solver.createSystem(dd);
  RealArr1D b("b", dd.ownedRowCount());
  RealArr1D x("x", dd.ownedRowCount());
  assembleLaplacian(dd, *sys, b);
  solver.solve(*sys, b, x);

  REQUIRE(solver.getIterationCount() >= 0);
  REQUIRE(solver.getResidualNorm() >= 0.0);
  REQUIRE(l2ErrorVsExact(dd, x) < 1e-8);  // rank-equivalence: matches exact on 1 and 2 ranks
}

TEST_CASE("pattern reuse: zeroEntries + reassemble + solve twice") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(8, 6, pnx, pny);
  Decomp2D dd(mc);

  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "gamg";
  cfg.rtol = 1e-12;
  PetscLinearSolver solver(cfg);
  auto sys = solver.createSystem(dd);
  RealArr1D b("b", dd.ownedRowCount());
  RealArr1D x1("x1", dd.ownedRowCount());
  RealArr1D x2("x2", dd.ownedRowCount());

  assembleLaplacian(dd, *sys, b);
  solver.solve(*sys, b, x1);
  // Reuse the same system/pattern without reallocating.
  assembleLaplacian(dd, *sys, b);  // calls zeroEntries internally
  solver.solve(*sys, b, x2);

  auto h1 = Kokkos::create_mirror_view(x1);
  auto h2 = Kokkos::create_mirror_view(x2);
  Kokkos::deep_copy(h1, x1);
  Kokkos::deep_copy(h2, x2);
  for (int L = 0; L < dd.ownedRowCount(); ++L)
    REQUIRE(h1(static_cast<size_t>(L)) == Approx(h2(static_cast<size_t>(L))).margin(1e-12));
  REQUIRE(l2ErrorVsExact(dd, x2) < 1e-8);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rc = frehg2test::runAll();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();   // does not finalize MPI (we called MPI_Init ourselves)
  MPI_Finalize();
  Kokkos::finalize();
  return grc;
}
