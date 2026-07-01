// P8.3.3: implicit diffusion in a 1D subsurface column.
//   - A step initial profile diffuses toward the analytical erfc solution.
//   - The Neumann operator conserves the column mass (sum of concentration).
//   - D <= 0 is a no-op.
// The implicit matrix is assembled and solved through the backend-agnostic LinearSolver /
// SparseSystem seam (PETSc is initialized here only to provide the backend).
#include <petscksp.h>

#include <cmath>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "frehg2/linear/LinearSolverFactory.hpp"
#include "frehg2/linear/SolverConfig.hpp"
#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "solute/Diffusion.hpp"

using namespace frehg2;

TEST_CASE("diffusion: 1D column matches the analytical erfc solution") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // 1D column reference is single-rank

  const int nz = 200;
  const real dz = 0.01;
  Grid grid(1, 1, nz, 1.0, 1.0, dz);
  MpiComm mc(1, 1, 1, 1);
  Decomp3D dd(mc, nz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-12;
  auto solver = makeLinearSolver("petsc", cfg);
  DiffusionSolver ds(*solver, dd, grid);

  // Step IC: C = 1 for k < nz/2 (z < z0), 0 below. z increases downward.
  RealArr1DHost c("c", grid.nCell());
  const int kstep = nz / 2;
  const real z0 = kstep * dz;
  for (int k = 0; k < nz; ++k)
    c(static_cast<size_t>(grid.getIndex(0, 0, k))) = (k < kstep) ? 1.0 : 0.0;

  const real D = 1.0e-5;
  const real dt = 0.25;
  const int nsteps = 250;
  const real t = dt * nsteps;  // = 62.5 s
  for (int s = 0; s < nsteps; ++s) ds.solve(c, D, dt);

  // Analytical: C(z,t) = 0.5*erfc((z - z0)/(2 sqrt(D t))).
  const real denom = 2.0 * std::sqrt(D * t);
  real max_err = 0.0;
  for (int k = 0; k < nz; ++k) {
    const real z = (k + 0.5) * dz;
    const real analytic = 0.5 * std::erfc((z - z0) / denom);
    const real model = c(static_cast<size_t>(grid.getIndex(0, 0, k)));
    max_err = std::max(max_err, std::fabs(model - analytic));
  }
  std::fprintf(stderr, "  diffusion erfc max abs error = %.3e\n", max_err);
  REQUIRE(max_err < 2.0e-2);
}

TEST_CASE("diffusion: Neumann operator conserves column mass") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nz = 120;
  const real dz = 0.01;
  Grid grid(1, 1, nz, 1.0, 1.0, dz);
  MpiComm mc(1, 1, 1, 1);
  Decomp3D dd(mc, nz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  auto solver = makeLinearSolver("petsc", cfg);
  DiffusionSolver ds(*solver, dd, grid);

  RealArr1DHost c("c", grid.nCell());
  for (int k = 0; k < nz; ++k)
    c(static_cast<size_t>(grid.getIndex(0, 0, k))) = (k < nz / 2) ? 1.0 : 0.0;
  real m0 = 0.0;
  for (int k = 0; k < nz; ++k) m0 += c(static_cast<size_t>(grid.getIndex(0, 0, k)));

  for (int s = 0; s < 50; ++s) ds.solve(c, 5.0e-5, 0.5);
  real m1 = 0.0;
  for (int k = 0; k < nz; ++k) m1 += c(static_cast<size_t>(grid.getIndex(0, 0, k)));
  REQUIRE(std::fabs(m1 - m0) < 1e-9 * std::fabs(m0));
}

TEST_CASE("diffusion: D <= 0 is a no-op") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nz = 10;
  Grid grid(1, 1, nz, 1.0, 1.0, 0.01);
  MpiComm mc(1, 1, 1, 1);
  Decomp3D dd(mc, nz);
  SolverConfig cfg;
  auto solver = makeLinearSolver("petsc", cfg);
  DiffusionSolver ds(*solver, dd, grid);

  RealArr1DHost c("c", grid.nCell());
  for (int k = 0; k < nz; ++k) c(static_cast<size_t>(grid.getIndex(0, 0, k))) = k + 1.0;
  ds.solve(c, 0.0, 1.0);
  for (int k = 0; k < nz; ++k)
    REQUIRE(c(static_cast<size_t>(grid.getIndex(0, 0, k))) == Approx(k + 1.0).margin(0.0));
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
