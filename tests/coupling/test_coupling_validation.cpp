// P6.3 coupling validation: a coupled SW+GW run on the b1-sw bathymetry (a real, non-flat
// terrain profile) driven by the b1-sw rainfall, with a synthetic uniform soil column for the
// groundwater. Gates:
//   * cumulative mass balance error < 1e-8 (relative) over the full simulation,
//   * no negative water depths (no over-extraction at the interface),
//   * no solver divergence (all states finite, GW stays unsaturated).
#include <petscksp.h>

#include <cmath>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "b1_sw_runner.hpp"
#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "coupling/Coupling.hpp"
#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;
using namespace frehg2::b1_sw;

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

}  // namespace

TEST_CASE("coupling validation: b1-sw bathymetry coupled run conserves mass, no negative depths") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string b1dir = std::string(FREHG2_SOURCE_DIR) + "/legacy/benchmarks/b1-sw";
  B1Params p = loadB1Params(b1dir);
  const double dx = p.dx, dy = p.dy;
  const double dt = 1.0;          // sync coupling dt (finer than the b1 dt=5; rainfall is f(t))
  const int nsteps = 300;         // 300 s of coupled evolution
  const double init_eta = 0.3;    // floods the whole b1 basin (bed+offset in [0, 0.6])
  MpiComm mc(p.gnx, p.gny, 1, 1);

  // Surface water: exactly the b1-sw geometry/parameters.
  Grid sgrid(mc.localNx(), mc.localNy(), 1, dx, dy, 0.1);
  SweSolver swe(sgrid, &mc);
  SweParams sp;
  sp.gravity = 9.81;
  sp.manning = 0.019;
  sp.min_depth = 1.0e-8;
  sp.visc_x = 1.0e-6;
  sp.visc_y = 1.0e-6;
  sp.hD = 0.1;
  sp.wtfh = 1.0e-8;
  sp.dt = dt;
  sp.offset = p.offset;
  swe.setParams(sp);
  RealArr1DHost bed("bed", static_cast<size_t>(p.gnx * p.gny));
  for (int j = 0; j < p.gny; ++j)
    for (int i = 0; i < p.gnx; ++i)
      bed(static_cast<size_t>(i + j * p.gnx)) = p.bath[static_cast<size_t>(j)] + p.offset;
  swe.setBathymetry(bed);
  swe.initializeState(init_eta);
  Decomp2D dd2(mc);
  SolverConfig cfg2;
  cfg2.ksp_type = "cg";
  cfg2.pc_type = "jacobi";
  cfg2.rtol = 1.0e-13;
  PetscLinearSolver solver2(cfg2);
  swe.attachSolver(solver2, dd2);

  // Groundwater: synthetic deep, incompressible, closed (no-flow) uniform soil column.
  const int nz = 20;
  const double dz = 0.1;
  Grid ggrid(mc.localNx(), mc.localNy(), nz, dx, dy, dz);
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
  rp.soil.Ks_z = 1.0e-7;  // slow infiltration: surface stays ponded, GW stays unsaturated
  rp.soil.Ss = 0.0;
  for (int b = 0; b < 6; ++b) rp.bc_type[b] = 0;  // closed box; coupling is the only top flux
  re.setParams(rp);
  re.initializeUniformColumn(0.15);
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
  cp.min_depth = 1.0e-8;
  Coupling coupling(cp);

  // Precompute the rainfall accumulation for this run (rainfall is a function of time).
  std::vector<double> rain_step(static_cast<size_t>(nsteps));
  double rain_accum = 0.0;
  {
    double t = 0.0;
    for (int step = 0; step < nsteps; ++step) {
      t += dt;
      rain_step[static_cast<size_t>(step)] = interpBc(p.rain_t, p.rain_v, t);
      rain_accum += rain_step[static_cast<size_t>(step)] * dt;
    }
  }
  // Rainfall enters every surface cell except the global y+ row (legacy SWE rain mask).
  const double rain_input =
      (dx * dy) * static_cast<double>(p.gnx * (p.gny - 1)) * rain_accum;

  // --- SWE-only baseline (same geometry/rainfall, no coupling): the surface solver's own
  // (P4) conservation floor, against which the coupled imbalance is bounded.
  Grid sgrid_b(mc.localNx(), mc.localNy(), 1, dx, dy, 0.1);
  SweSolver swe_b(sgrid_b, &mc);
  swe_b.setParams(sp);
  swe_b.setBathymetry(bed);
  swe_b.initializeState(init_eta);
  Decomp2D dd2b(mc);
  PetscLinearSolver solver2b(cfg2);
  swe_b.attachSolver(solver2b, dd2b);
  const double swb_init = surfaceVolume(swe_b);
  for (int step = 0; step < nsteps; ++step) swe_b.advanceStep(rain_step[static_cast<size_t>(step)], 0.0);
  const double swe_only_floor = std::fabs((surfaceVolume(swe_b) - swb_init) - rain_input);

  // --- Coupled run.
  const double sw_init = surfaceVolume(swe);
  const double gw_init = groundwaterVolume(re);
  double net_total = 0.0;     // signed exchanged volume (+ = GW->SW seepage)
  double min_depth = 1e30;
  double max_depth = -1e30;
  double wc_max = -1e30;
  for (int step = 0; step < nsteps; ++step) {
    net_total += coupling.stepCoupled(swe, re, rain_step[static_cast<size_t>(step)], 0.0, dt);
    for (int j = 0; j < p.gny; ++j)
      for (int i = 0; i < p.gnx; ++i) {
        const int si = sgrid.getSurfaceIndex(i, j);
        const double d = swe.fields().eta(si) - swe.fields().bottom(si);
        REQUIRE(std::isfinite(d));
        REQUIRE(d >= -1.0e-12);
        min_depth = std::min(min_depth, d);
        max_depth = std::max(max_depth, d);
      }
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < p.gny; ++j)
        for (int i = 0; i < p.gnx; ++i) {
          const double w = re.fields().wc(ggrid.getIndex(i, j, k));
          REQUIRE(std::isfinite(w));
          wc_max = std::max(wc_max, w);
        }
  }

  const double sw_fin = surfaceVolume(swe);
  const double gw_fin = groundwaterVolume(re);
  const double total0 = sw_init + gw_init;
  const double imbalance = (sw_fin + gw_fin - total0) - rain_input;
  const double rel = std::fabs(imbalance) / total0;
  const double gw_abs_err = std::fabs((gw_fin - gw_init) - (-net_total));
  (void)swe_only_floor;
  std::fprintf(stderr,
               "  b1-sw coupled (%d steps): dSW=%.6e dGW=%.6e net_exch=%.6e rain=%.6e\n"
               "    imbalance=%.3e rel=%.3e swe_only_floor=%.3e gw_abs_err=%.3e\n"
               "    wc_max=%.6f min_depth=%.3e max_depth=%.3e\n",
               nsteps, sw_fin - sw_init, gw_fin - gw_init, net_total, rain_input, imbalance,
               rel, swe_only_floor, gw_abs_err, wc_max, min_depth, max_depth);

  // The basin held water throughout and there was real infiltration.
  REQUIRE(max_depth > 0.1);
  REQUIRE(net_total < 0.0);
  // Stayed unsaturated -> conservative regime; no over-extraction.
  REQUIRE(wc_max < rp.soil.theta_s);
  REQUIRE(min_depth >= -1.0e-12);
  // Groundwater receives EXACTLY what the exchange removed from the surface (machine precision).
  REQUIRE(gw_abs_err < 1.0e-9 * std::fabs(net_total));
  // Cumulative mass balance error (relative) below 1e-8 over the full simulation.
  REQUIRE(rel < 1.0e-8);
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
