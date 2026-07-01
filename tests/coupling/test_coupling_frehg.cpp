// P6.1 coupling-mechanism gate: SW<->GW exchange-flux sign conventions, K_sat-limited
// infiltration, and full-grid (not single-column) coverage. The exchange physics is the
// z-back top-face form of legacy darcy_flux (subroutines.c:70-112) with the plan's total-head
// sign convention (P6.1): q>0 seepage (GW->SW), q<0 infiltration (SW->GW), q==0 at equal head.
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

namespace {
CouplingParams makeParams() {
  CouplingParams cp;
  cp.dx = 10.0;
  cp.dy = 10.0;
  cp.visc = 1.0;
  cp.min_depth = 1.0e-8;
  return cp;
}
}  // namespace

TEST_CASE("coupling: exchange-flux sign conventions (infiltration/seepage/zero)") {
  Coupling coupling(makeParams());
  const double Ksz = 1.0e-5, dz_top = 0.1, dept = 0.5;

  // h_surface > h_gw_top -> infiltration (SW -> GW), q < 0.
  const double q_inf = coupling.columnExchangeRate(/*h_gw_top=*/0.0, /*h_surface=*/1.0, dept,
                                                   Ksz, dz_top);
  REQUIRE(q_inf < 0.0);

  // h_gw_top > h_surface -> seepage (GW -> SW), q > 0.
  const double q_seep = coupling.columnExchangeRate(/*h_gw_top=*/1.0, /*h_surface=*/0.0, dept,
                                                    Ksz, dz_top);
  REQUIRE(q_seep > 0.0);

  // Equal heads -> exactly zero flux (no spurious gravity term).
  const double q_zero = coupling.columnExchangeRate(0.5, 0.5, dept, Ksz, dz_top);
  REQUIRE(q_zero == Approx(0.0).margin(1.0e-30));

  // Antisymmetry in the head gradient.
  REQUIRE(q_inf == Approx(-q_seep).margin(1.0e-30));

  // A second parameter set (different soil/geometry) keeps the same signs.
  const double q_inf2 = coupling.columnExchangeRate(0.2, 0.9, 0.3, 5.0e-6, 0.05);
  const double q_seep2 = coupling.columnExchangeRate(0.9, 0.2, 0.3, 5.0e-6, 0.05);
  REQUIRE(q_inf2 < 0.0);
  REQUIRE(q_seep2 > 0.0);
}

TEST_CASE("coupling: infiltration capacity scales with K_sat, no ponded water -> no infiltration") {
  const CouplingParams cp = makeParams();
  Coupling coupling(cp);
  const double dz_top = 0.1, dept = 0.5, dh = 1.0;  // h_gw_top - h_surface
  const double Az = cp.dx * cp.dy;
  const double delta = 0.5 * dz_top;

  // Exact closed form: q = Az * Ksz * visc * dh / delta.
  const double Ksz1 = 1.0e-5;
  const double q1 = coupling.columnExchangeRate(dh, 0.0, dept, Ksz1, dz_top);
  REQUIRE(q1 == Approx(Az * Ksz1 * cp.visc * dh / delta).epsilon(1.0e-12));

  // K_sat sets the conductance: doubling Ksz doubles the flux magnitude.
  const double q2 = coupling.columnExchangeRate(dh, 0.0, dept, 2.0 * Ksz1, dz_top);
  REQUIRE(q2 == Approx(2.0 * q1).epsilon(1.0e-12));

  // No ponded surface water: an infiltration gradient (h_surface > h_gw_top) draws nothing.
  const double q_dry_inf = coupling.columnExchangeRate(0.0, 1.0, /*dept=*/0.0, Ksz1, dz_top);
  REQUIRE(q_dry_inf == Approx(0.0).margin(1.0e-30));

  // ... but a dry surface can still receive seepage when the GW head is higher.
  const double q_dry_seep = coupling.columnExchangeRate(1.0, 0.0, /*dept=*/0.0, Ksz1, dz_top);
  REQUIRE(q_dry_seep > 0.0);
}

TEST_CASE("coupling: exchange computed for every cell of a 10x10 grid") {
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

  // Spatially varying surface elevation and top-GW pressure head so each column differs.
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      swe.fields().eta(sgrid.getSurfaceIndex(i, j)) = 0.3 + 0.05 * i + 0.02 * j;
      re.fields().h(ggrid.getIndex(i, j, 0)) = -0.5 + 0.03 * j - 0.01 * i;
    }

  RealArr1DHost q("q", static_cast<size_t>(nx * ny));
  CouplingParams cp = makeParams();
  cp.dx = dx;
  cp.dy = dy;
  Coupling coupling(cp);
  coupling.computeExchangeRates(swe, re, q);

  // Every column must carry the value predicted by the per-column closed form, and the
  // field must be genuinely spatially varying (not one hard-coded column).
  const SoilParams& soil = re.params().soil;
  double qmin = 1e30, qmax = -1e30;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int si = sgrid.getSurfaceIndex(i, j);
      const int ti = ggrid.getIndex(i, j, 0);
      const double dept = std::max(swe.fields().eta(si) - swe.fields().bottom(si), 0.0);
      const double dz_top = re.fields().dz3d(ti);
      const double h_gw = re.fields().h(ti);  // top-cell pressure head
      const double expect =
          coupling.columnExchangeRate(h_gw, dept, dept, soil.Ks_z, dz_top);
      const double got = q(static_cast<size_t>(i + j * nx));
      REQUIRE(got == Approx(expect).epsilon(1.0e-13).margin(1.0e-300));
      qmin = std::min(qmin, got);
      qmax = std::max(qmax, got);
    }
  // Field genuinely varies across the grid.
  REQUIRE(qmax - qmin > 1.0e-9);
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
