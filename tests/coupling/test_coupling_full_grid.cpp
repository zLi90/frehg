// P6.1 full-grid coupling gate: on a 10x10 grid with spatially varying surface and GW state,
// each column's exchange flux must depend ONLY on that column's state. Perturbing one column
// changes that column's flux and leaves every other column bit-identical (no cross-column
// coupling in the exchange computation).
#include <cmath>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "coupling/Coupling.hpp"
#include "frehg2_test.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

TEST_CASE("coupling full grid: per-column independence (no cross-column coupling)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 10, ny = 10, nz = 5;
  const double dx = 10.0, dy = 10.0, dz = 0.1;
  MpiComm mc(nx, ny, 1, 1);

  Grid sgrid(nx, ny, 1, dx, dy, 1.0);
  SweSolver swe(sgrid, &mc);
  SweParams sp;
  sp.dt = 1.0;
  swe.setParams(sp);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(0.5);

  Grid ggrid(nx, ny, nz, dx, dy, dz);
  ReSolver re(ggrid, &mc);
  ReParams rp;
  rp.dx = dx;
  rp.dy = dy;
  rp.dz = dz;
  re.setParams(rp);
  re.initializeUniformColumn(0.1);

  // Spatially varying surface elevation and top-GW pressure head.
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      swe.fields().eta(sgrid.getSurfaceIndex(i, j)) =
          0.4 + 0.03 * i - 0.02 * j + 0.001 * i * j;
      re.fields().h(ggrid.getIndex(i, j, 0)) = -0.4 + 0.02 * j - 0.015 * i;
    }

  CouplingParams cp;
  cp.dx = dx;
  cp.dy = dy;
  cp.visc = 1.0;
  cp.min_depth = 1.0e-8;
  Coupling coupling(cp);

  RealArr1DHost q0("q0", static_cast<size_t>(nx * ny));
  coupling.computeExchangeRates(swe, re, q0);

  // Perturb a single interior column's surface elevation and top GW head.
  const int pi = 4, pj = 7;
  swe.fields().eta(sgrid.getSurfaceIndex(pi, pj)) += 0.25;
  re.fields().h(ggrid.getIndex(pi, pj, 0)) += 0.17;

  RealArr1DHost q1("q1", static_cast<size_t>(nx * ny));
  coupling.computeExchangeRates(swe, re, q1);

  // Only the perturbed column changes; all others are bit-identical.
  int changed = 0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const size_t oi = static_cast<size_t>(i + j * nx);
      if (i == pi && j == pj) {
        REQUIRE(std::fabs(q1(oi) - q0(oi)) > 1.0e-9);
        ++changed;
      } else {
        REQUIRE(q1(oi) == q0(oi));  // exact bitwise equality off the perturbed column
      }
    }
  REQUIRE(changed == 1);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  int rc = frehg2test::runAll();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  Kokkos::finalize();
  return grc;
}
