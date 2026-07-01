// P8.3.5: the standalone SoluteStepper driven by a synthetic flow field (no Orchestrator).
//   - solute.enabled == false is a no-op.
//   - Advection on a closed domain conserves mass to machine precision.
//   - Decay drives the total mass as exp(-k t).
//   - The CFL safety check refuses the step and leaves the state unchanged.
//   - The full source->advect->diffuse->record pipeline (with implicit diffusion attached)
//     conserves mass.
#include <petscksp.h>

#include <cmath>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "core/Grid.hpp"
#include "core/GwState.hpp"
#include "core/MpiComm.hpp"
#include "core/State.hpp"
#include "frehg2/linear/LinearSolverFactory.hpp"
#include "frehg2/linear/SolverConfig.hpp"
#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "solute/SoluteFlow.hpp"
#include "solute/SoluteParams.hpp"
#include "solute/SoluteStepper.hpp"

using namespace frehg2;

namespace {

real gaussian(real x, real x0, real sd) {
  const real z = (x - x0) / sd;
  return std::exp(-0.5 * z * z);
}

void setSurfaceGaussian(State& sw, const Grid& g, real x0, real sd) {
  auto h = Kokkos::create_mirror_view(sw.conc);
  Kokkos::deep_copy(h, sw.conc);
  for (int j = 0; j < g.ny(); ++j)
    for (int i = 0; i < g.nx(); ++i)
      h(static_cast<size_t>(g.getSurfaceIndex(i, j))) = gaussian(i, x0, sd);
  Kokkos::deep_copy(sw.conc, h);
}

real surfaceMass(const State& sw, const Grid& g) {
  auto h = Kokkos::create_mirror_view(sw.conc);
  Kokkos::deep_copy(h, sw.conc);
  real s = 0.0;
  for (int j = 0; j < g.ny(); ++j)
    for (int i = 0; i < g.nx(); ++i)
      s += h(static_cast<size_t>(g.getSurfaceIndex(i, j)));
  return s;
}

SoluteFlow makeUniformSurfaceFlow(const Grid& g, real u0) {
  SoluteFlow flow;
  flow.u = RealArr1DHost("u", g.nSurfaceStorageCell());
  flow.v = RealArr1DHost("v", g.nSurfaceStorageCell());
  flow.depth = RealArr1DHost("d", g.nSurfaceStorageCell());
  for (size_t i = 0; i < flow.u.extent(0); ++i) {
    flow.u(i) = u0;
    flow.v(i) = 0.0;
    flow.depth(i) = 1.0;
  }
  return flow;
}

}  // namespace

TEST_CASE("stepper: disabled solute is a no-op") {
  Grid grid(16, 1, 1, 1.0, 1.0, 1.0);
  State sw(grid);
  GwState gw(grid);
  setSurfaceGaussian(sw, grid, 8.0, 2.0);
  const real m0 = surfaceMass(sw, grid);
  SoluteParams p;  // enabled == false
  SoluteStepper stepper(grid, p);
  SoluteFlow flow = makeUniformSurfaceFlow(grid, 1.0);
  StepResult r = stepper.step(sw, gw, flow, 1.0);
  REQUIRE(r.ok);
  REQUIRE(surfaceMass(sw, grid) == Approx(m0).margin(0.0));
}

TEST_CASE("stepper: advection conserves mass on a closed domain") {
  Grid grid(30, 1, 1, 1.0, 1.0, 1.0);
  State sw(grid);
  GwState gw(grid);
  setSurfaceGaussian(sw, grid, 8.0, 2.0);
  const real m0 = surfaceMass(sw, grid);
  SoluteParams p;
  p.enabled = true;
  p.diffusion_scheme = "none";
  SoluteStepper stepper(grid, p);
  SoluteFlow flow = makeUniformSurfaceFlow(grid, 0.6);
  for (int s = 0; s < 12; ++s) {
    StepResult r = stepper.step(sw, gw, flow, 1.0);
    REQUIRE(r.ok);
    REQUIRE(r.max_cfl == Approx(0.6).margin(1e-12));
  }
  REQUIRE(std::fabs(surfaceMass(sw, grid) - m0) < 1e-12 * std::fabs(m0));
}

TEST_CASE("stepper: decay drives total mass as exp(-k t)") {
  Grid grid(30, 1, 1, 1.0, 1.0, 1.0);
  State sw(grid);
  GwState gw(grid);
  setSurfaceGaussian(sw, grid, 10.0, 2.5);
  const real m0 = surfaceMass(sw, grid);
  SoluteParams p;
  p.enabled = true;
  p.diffusion_scheme = "none";
  p.k_decay = 0.05;
  SoluteStepper stepper(grid, p);
  SoluteFlow flow = makeUniformSurfaceFlow(grid, 0.5);
  const real dt = 1.0;
  const int n = 15;
  for (int s = 0; s < n; ++s) {
    StepResult r = stepper.step(sw, gw, flow, dt);
    REQUIRE(r.ok);
  }
  const real expected = m0 * std::exp(-p.k_decay * dt * n);
  REQUIRE(surfaceMass(sw, grid) == Approx(expected).epsilon(1e-10));
}

TEST_CASE("stepper: CFL refusal leaves the state unchanged") {
  Grid grid(16, 1, 1, 1.0, 1.0, 1.0);
  State sw(grid);
  GwState gw(grid);
  setSurfaceGaussian(sw, grid, 8.0, 2.0);
  const real m0 = surfaceMass(sw, grid);
  SoluteParams p;
  p.enabled = true;
  p.diffusion_scheme = "none";
  SoluteStepper stepper(grid, p);
  SoluteFlow flow = makeUniformSurfaceFlow(grid, 1.0);
  StepResult r = stepper.step(sw, gw, flow, 1.5);  // Courant 1.5 > 1
  REQUIRE_FALSE(r.ok);
  REQUIRE(r.max_cfl > p.cfl_max);
  REQUIRE(surfaceMass(sw, grid) == Approx(m0).margin(0.0));
}

TEST_CASE("stepper: full advect+diffuse pipeline conserves mass") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  Grid grid(24, 1, 1, 1.0, 1.0, 1.0);
  State sw(grid);
  GwState gw(grid);
  setSurfaceGaussian(sw, grid, 8.0, 2.0);
  const real m0 = surfaceMass(sw, grid);

  SoluteParams p;
  p.enabled = true;
  p.diffusion_scheme = "implicit";
  p.D = 1.0e-3;
  SoluteStepper stepper(grid, p);

  MpiComm mc(grid.nx(), grid.ny(), 1, 1);
  Decomp2D dd(mc);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  auto solver = makeLinearSolver("petsc", cfg);
  stepper.attachSurfaceDiffusion(*solver, dd);

  SoluteFlow flow = makeUniformSurfaceFlow(grid, 0.5);
  for (int s = 0; s < 10; ++s) {
    StepResult r = stepper.step(sw, gw, flow, 1.0);
    REQUIRE(r.ok);
  }
  REQUIRE(std::fabs(surfaceMass(sw, grid) - m0) < 1e-8 * std::fabs(m0));
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
