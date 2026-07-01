// P4.0 b0-lake well-balanced ("lake at rest") test, serial one rank.
//
// A flat free surface (eta=1.0) over a Gaussian bed bump must remain exactly at rest under
// the semi-implicit SWE scheme: no spurious velocity, no surface drift. This exercises the
// momentum/RHS/matrix/solve/velocity chain (fully wet, no dry-cell hack) through the
// SparseSystem/LinearSolver seam.
#include <petscksp.h>

#include <cmath>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "core/MpiComm.hpp"
#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

TEST_CASE("b0-lake: flat surface over a bump stays at rest (well-balanced)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 20, ny = 20;
  const double dx = 1.0, dy = 1.0, dt = 1.0;
  const double eta0 = 1.0;
  const int nsteps = 200;

  MpiComm mc(nx, ny, 1, 1);
  Grid grid(mc.localNx(), mc.localNy(), 1, dx, dy, 1.0);
  SweSolver swe(grid, &mc);
  SweParams p;
  p.gravity = 9.81;
  p.manning = 0.02;
  p.min_depth = 1.0e-8;
  p.hD = 0.1;
  p.wtfh = 1.0e-8;
  p.dt = dt;
  p.offset = 0.0;
  swe.setParams(p);

  // Gaussian bed bump z = 0.5 exp(-((x-10)^2 + (y-10)^2)/10), cell-centered.
  RealArr1DHost bed("bed", static_cast<size_t>(nx * ny));
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const double x = (i + 0.5) * dx, y = (j + 0.5) * dy;
      bed(static_cast<size_t>(i + j * nx)) =
          0.5 * std::exp(-((x - 10.0) * (x - 10.0) + (y - 10.0) * (y - 10.0)) / 10.0);
    }
  swe.setBathymetry(bed);
  swe.initializeState(eta0);

  Decomp2D dd(mc);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1e-13;
  PetscLinearSolver solver(cfg);
  swe.attachSolver(solver, dd);

  for (int step = 0; step < nsteps; ++step) swe.advanceStep(0.0, 0.0);

  double max_eta_dev = 0.0, max_vel = 0.0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = grid.getSurfaceIndex(i, j);
      max_eta_dev = std::max(max_eta_dev, std::fabs(swe.fields().eta(c) - eta0));
      max_vel = std::max(max_vel, std::fabs(swe.fields().uu(c)));
      max_vel = std::max(max_vel, std::fabs(swe.fields().vv(c)));
    }
  std::fprintf(stderr, "  b0-lake after %d steps: max|eta-1| = %.3e, max|vel| = %.3e\n",
               nsteps, max_eta_dev, max_vel);
  REQUIRE(max_eta_dev < 1.0e-8);
  REQUIRE(max_vel < 1.0e-8);
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
