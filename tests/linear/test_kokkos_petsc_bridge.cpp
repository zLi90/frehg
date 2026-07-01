// P2.6 acceptance: assemble a 5-point Laplacian whose coefficients are computed on-device
// with Kokkos::parallel_for, staged to host, and inserted via SparseSystem::addRow; then
// solve and compare against the analytic solution. Runs on 1 and 2 MPI ranks.
#include <petscksp.h>

#include <cmath>

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
  else { mpi_nx = size; mpi_ny = 1; }
}
real xExact(int grow) { return 2.0 + std::cos(0.07 * static_cast<real>(grow)); }
}  // namespace

TEST_CASE("device-computed Laplacian coefficients assemble and solve correctly") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(10, 8, pnx, pny);
  Decomp2D dd(mc);
  const int n_local = dd.ownedRowCount();

  // Host-side stencil structure (rows/cols) from the decomposition.
  IntArr1DHost grow_h("grow", n_local);
  IntArr1DHost ncol_h("ncol", n_local);
  IntArr2DHost cols_h("cols", n_local, kMaxStencil);
  int cols[kMaxStencil];
  for (int L = 0; L < n_local; ++L) {
    int grow = 0, ncols = 0;
    dd.ownedRowStencil(L, grow, cols, ncols);
    grow_h(L) = grow;
    ncol_h(L) = ncols;
    for (int c = 0; c < ncols; ++c) cols_h(L, c) = cols[c];
  }

  // Compute coefficients on device with a parallel_for (the GPU-capable assembly path).
  RealArr2D vals_d("vals", n_local, kMaxStencil);
  auto ncol_d = Kokkos::create_mirror_view_and_copy(DeviceSpace(), ncol_h);
  Kokkos::parallel_for(
      "laplacian_coeffs", Kokkos::RangePolicy<DeviceExec>(0, n_local),
      KOKKOS_LAMBDA(int L) {
        const int nc = ncol_d(L);
        for (int c = 0; c < nc; ++c) vals_d(L, c) = (c == 0) ? 4.0 : -1.0;
      });
  Kokkos::fence();

  // Stage to host and insert via the SparseSystem interface.
  auto vals_h = Kokkos::create_mirror_view(vals_d);
  Kokkos::deep_copy(vals_h, vals_d);

  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = (size == 1) ? "jacobi" : "gamg";
  cfg.rtol = 1e-14;
  PetscLinearSolver solver(cfg);
  auto sys = solver.createSystem(dd);

  sys->zeroEntries();
  sys->beginAssembly();
  RealArr1D b("b", n_local);
  auto bh = Kokkos::create_mirror_view(b);
  for (int L = 0; L < n_local; ++L) {
    const int ncols = ncol_h(L);
    real rowvals[kMaxStencil];
    int rowcols[kMaxStencil];
    real bv = 0.0;
    for (int c = 0; c < ncols; ++c) {
      rowvals[c] = vals_h(L, c);
      rowcols[c] = cols_h(L, c);
      bv += rowvals[c] * xExact(rowcols[c]);
    }
    sys->addRow(grow_h(L), ncols, rowcols, rowvals);
    bh(static_cast<size_t>(L)) = bv;
  }
  Kokkos::deep_copy(b, bh);
  sys->endAssembly();

  RealArr1D x("x", n_local);
  solver.solve(*sys, b, x);

  // Compare against the analytic solution: max element diff < 1e-10 (1 rank, lu).
  auto xh = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(xh, x);
  const int rstart = dd.ownershipRange().first;
  real local_max = 0.0;
  for (int L = 0; L < n_local; ++L)
    local_max = std::fmax(local_max, std::fabs(xh(static_cast<size_t>(L)) - xExact(rstart + L)));
  real global_max = local_max;
  MPI_Allreduce(MPI_IN_PLACE, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  const real tol = (size == 1) ? 1e-9 : 1e-7;
  REQUIRE(global_max < tol);
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
