// P6.2 synchronous coupled time-stepping gate.
//   * The exchange is exactly mass-conservative: the volume added to one domain equals the
//     volume removed from the other (machine precision, independent of solver tolerances).
//   * A coupled SW+GW run on a full 10x10x5 grid is stable for 100 steps (finite depths >=0,
//     water content within [theta_r, theta_s]), conserves total water to the rainfall input,
//     and the exchange direction follows the head gradient.
#include <petscksp.h>

#include <cmath>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "coupling/Coupling.hpp"
#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

namespace {

double surfaceVolume(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const double Az = g.dx() * g.dy();
  double v = 0.0;
  for (int j = 0; j < g.ny(); ++j)
    for (int i = 0; i < g.nx(); ++i) {
      const int si = g.getSurfaceIndex(i, j);
      v += std::max(swe.fields().eta(si) - swe.fields().bottom(si), 0.0) * Az;
    }
  return v;
}

double groundwaterVolume(const ReSolver& re) {
  const Grid& g = re.grid();
  const double Az = g.dx() * g.dy();
  double v = 0.0;
  for (int k = 0; k < g.nz(); ++k)
    for (int j = 0; j < g.ny(); ++j)
      for (int i = 0; i < g.nx(); ++i) {
        const int c = g.getIndex(i, j, k);
        v += re.fields().wc(c) * Az * re.fields().dz3d(c);
      }
  return v;
}

double maxWaterContent(const ReSolver& re) {
  const Grid& g = re.grid();
  double m = -1e30;
  for (int k = 0; k < g.nz(); ++k)
    for (int j = 0; j < g.ny(); ++j)
      for (int i = 0; i < g.nx(); ++i)
        m = std::max(m, re.fields().wc(g.getIndex(i, j, k)));
  return m;
}

}  // namespace

TEST_CASE("coupled: exchange is exactly mass-conservative (SW gain == GW loss)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 10, ny = 10, nz = 5;
  const double dx = 10.0, dy = 10.0, dz = 0.1, dt = 1.0;
  MpiComm mc(nx, ny, 1, 1);

  Grid sgrid(nx, ny, 1, dx, dy, 1.0);
  SweSolver swe(sgrid, &mc);
  SweParams sp;
  sp.dt = dt;
  swe.setParams(sp);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(0.5);

  Grid ggrid(nx, ny, nz, dx, dy, dz);
  ReSolver re(ggrid, &mc);
  ReParams rp;
  rp.dx = dx;
  rp.dy = dy;
  rp.dz = dz;
  rp.soil.theta_s = 0.33;
  rp.soil.theta_r = 0.0;
  re.setParams(rp);
  re.initializeUniformColumn(0.2);

  CouplingParams cp;
  cp.dx = dx;
  cp.dy = dy;
  Coupling coupling(cp);

  // A net-nonzero exchange field (varying-magnitude infiltration) so dSW is genuinely
  // nonzero, small enough that neither donor runs dry.
  RealArr1DHost q("q", static_cast<size_t>(nx * ny));
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      q(static_cast<size_t>(i + j * nx)) = -1.0e-4 * (1.0 + 0.5 * ((i + j) % 3));

  double q_total = 0.0;  // analytic transferred volume = sum(q*dt)
  for (int idx = 0; idx < nx * ny; ++idx) q_total += q(static_cast<size_t>(idx)) * dt;

  const double sw0 = surfaceVolume(swe);
  const double gw0 = groundwaterVolume(re);
  coupling.applyExchangeToSurface(q, swe, dt);
  coupling.applyExchangeToGroundwater(q, re, dt);
  const double dSW = surfaceVolume(swe) - sw0;
  const double dGW = groundwaterVolume(re) - gw0;
  std::fprintf(stderr, "  exchange: dSW=%.6e dGW=%.6e  sum=%.3e  q_total=%.6e\n", dSW, dGW,
               dSW + dGW, q_total);
  // Meaningful (nonzero) transfer; the surface gains and groundwater loses exactly q_total.
  // (Tolerance covers volume-differencing cancellation against the ~0.5 m baseline depth.)
  REQUIRE(std::fabs(q_total) > 1.0e-3);
  REQUIRE(std::fabs(dSW - q_total) < 1.0e-9);
  REQUIRE(std::fabs(dGW + q_total) < 1.0e-9);
  REQUIRE(std::fabs(dSW + dGW) < 1.0e-9);
}

TEST_CASE("coupled: SW+GW run is stable and conserves mass over 100 steps") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nx = 10, ny = 10, nz = 5;
  const double dx = 10.0, dy = 10.0, dz = 0.1;
  const double dt = 0.5;
  const int nsteps = 100;
  const double rain = 1.0e-6;  // m/s, applied to the surface only
  MpiComm mc(nx, ny, 1, 1);

  // Surface: closed basin, flat bed, shallow pond.
  Grid sgrid(nx, ny, 1, dx, dy, 1.0);
  SweSolver swe(sgrid, &mc);
  SweParams sp;
  sp.gravity = 9.81;
  sp.manning = 0.0;
  sp.min_depth = 1.0e-10;
  sp.dt = dt;
  swe.setParams(sp);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(0.1);
  Decomp2D dd2(mc);
  SolverConfig cfg2;
  cfg2.ksp_type = "cg";
  cfg2.pc_type = "jacobi";
  cfg2.rtol = 1.0e-14;
  PetscLinearSolver solver2(cfg2);
  swe.attachSolver(solver2, dd2);

  // Groundwater: closed box (all no-flow), incompressible (Ss=0 -> exact wc conservation),
  // fixed dt to stay synchronous with the surface, unsaturated throughout.
  Grid ggrid(nx, ny, nz, dx, dy, dz);
  ReSolver re(ggrid, &mc);
  ReParams rp;
  rp.dx = dx;
  rp.dy = dy;
  rp.dz = dz;
  rp.dt = dt;
  rp.adaptive_dt = false;
  rp.use_corrector = true;
  rp.soil.theta_s = 0.33;
  rp.soil.theta_r = 0.0;
  rp.soil.Ks_z = 2.89e-6;
  rp.soil.Ss = 0.0;
  for (int b = 0; b < 6; ++b) rp.bc_type[b] = 0;  // all no-flow; coupling is the only top flux
  re.setParams(rp);
  re.initializeUniformColumn(0.2);
  Decomp3D dd3(mc, nz);
  SolverConfig cfg3;
  cfg3.ksp_type = "cg";
  cfg3.pc_type = "sor";
  cfg3.rtol = 1.0e-12;
  PetscLinearSolver solver3(cfg3);
  re.attachSolver(solver3, dd3);

  CouplingParams cp;
  cp.dx = dx;
  cp.dy = dy;
  cp.min_depth = 1.0e-10;
  Coupling coupling(cp);

  // Rainfall enters every surface cell except the global y+ row (legacy SWE rain mask).
  const double rain_input =
      rain * dt * (dx * dy) * static_cast<double>(nx * (ny - 1)) * nsteps;

  // --- SWE-only baseline: isolates the surface solver's own (P4) conservation floor.
  Grid sgrid_b(nx, ny, 1, dx, dy, 1.0);
  SweSolver swe_b(sgrid_b, &mc);
  swe_b.setParams(sp);
  swe_b.setBathymetryConstant(0.0);
  swe_b.initializeState(0.1);
  Decomp2D dd2b(mc);
  PetscLinearSolver solver2b(cfg2);
  swe_b.attachSolver(solver2b, dd2b);
  const double swb_init = surfaceVolume(swe_b);
  for (int step = 0; step < nsteps; ++step) swe_b.advanceStep(rain, 0.0);
  const double swe_only_floor = std::fabs((surfaceVolume(swe_b) - swb_init) - rain_input);

  // --- Coupled run.
  const double sw_init = surfaceVolume(swe);
  const double gw_init = groundwaterVolume(re);

  // Initial head gradient: ponded surface over unsaturated soil -> net infiltration (q<0).
  const double net_first = coupling.stepCoupled(swe, re, rain, 0.0, dt);
  REQUIRE(net_first < 0.0);

  double net_total = net_first;
  double min_depth_seen = 1e30;
  for (int step = 1; step < nsteps; ++step) {
    net_total += coupling.stepCoupled(swe, re, rain, 0.0, dt);
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const double d = swe.fields().eta(sgrid.getSurfaceIndex(i, j)) -
                         swe.fields().bottom(sgrid.getSurfaceIndex(i, j));
        REQUIRE(std::isfinite(d));
        REQUIRE(d >= -1.0e-12);
        min_depth_seen = std::min(min_depth_seen, d);
      }
  }

  const double sw_fin = surfaceVolume(swe);
  const double gw_fin = groundwaterVolume(re);
  const double wc_max = maxWaterContent(re);

  // Stayed strictly unsaturated (so the non-conservative reallocation branch never fired).
  REQUIRE(wc_max < rp.soil.theta_s);

  const double total0 = sw_init + gw_init;
  const double imbalance = (sw_fin + gw_fin - total0) - rain_input;
  const double rel = std::fabs(imbalance) / total0;
  // GW absorbs EXACTLY the exchanged volume (-net_total): coupling+GW is conservative to
  // machine precision, independent of any solver tolerance.
  const double gw_abs_err = std::fabs((gw_fin - gw_init) - (-net_total));
  std::fprintf(stderr,
               "  coupled 100 steps: dSW=%.6e dGW=%.6e net_exch=%.6e rain=%.6e\n"
               "    imbalance=%.3e rel=%.3e  swe_only_floor=%.3e gw_abs_err=%.3e\n"
               "    wc_max=%.6f min_depth=%.3e\n",
               sw_fin - sw_init, gw_fin - gw_init, net_total, rain_input, imbalance, rel,
               swe_only_floor, gw_abs_err, wc_max, min_depth_seen);

  // Sustained infiltration (exchange direction follows the head gradient).
  REQUIRE(net_total < 0.0);
  // Groundwater receives exactly what the exchange removed from the surface (machine precision).
  REQUIRE(gw_abs_err < 1.0e-9 * std::fabs(net_total));
  // System mass balance: the ONLY non-exact contribution is the surface solver, so the total
  // imbalance must not exceed the SWE-only baseline floor (the coupling adds no mass loss).
  REQUIRE(std::fabs(imbalance) <= swe_only_floor + 1.0e-9);
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
