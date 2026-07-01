// P7.1 acceptance: the Orchestrator (the single production path) reproduces the P4 (b1-sw)
// and P5 (b2-gw) direct-call solver paths bit-for-bit. We hand-wire each direct run exactly
// as the solvers were exercised in P4/P5 and assert the unified-driver result matches to
// L2 < 1e-10 (same executable/backend).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <fstream>
#include <mpi.h>
#include <string>
#include <vector>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;
const std::string kB1Dir = std::string(FREHG2_SOURCE_DIR) + "/benchmarks/b1-sw";

std::vector<double> readDoubles(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<double> v;
  double x;
  while (in >> x) v.push_back(x);
  return v;
}

double interpBc(const std::vector<double>& t, const std::vector<double>& val, double tc) {
  const int n = static_cast<int>(t.size());
  if (tc == 0.0) return val[0];
  int ind = 1;
  if (ind < n && tc > t[ind]) {
    ind += 1;
    while (ind < n && tc > t[ind]) ind += 1;
  }
  if (ind >= n) return val[static_cast<size_t>(n - 1)];
  return val[static_cast<size_t>(ind - 1)] +
         (val[static_cast<size_t>(ind)] - val[static_cast<size_t>(ind - 1)]) *
             (tc - t[static_cast<size_t>(ind - 1)]) /
             (t[static_cast<size_t>(ind)] - t[static_cast<size_t>(ind - 1)]);
}

RealArr1DHost ownedDepth(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const int nx = g.nx(), ny = g.ny();
  RealArr1DHost d("d", static_cast<size_t>(nx * ny));
  const auto& f = swe.fields();
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = g.getSurfaceIndex(i, j);
      d(static_cast<size_t>(i + j * nx)) = std::max(f.eta(c) - f.bottom(c), 0.0);
    }
  return d;
}

RealArr1DHost wcColumn(const ReSolver& re) {
  const Grid& g = re.grid();
  const int nz = g.nz();
  RealArr1DHost w("w", static_cast<size_t>(nz));
  for (int k = 0; k < nz; ++k) w(static_cast<size_t>(k)) = re.fields().wc(g.getIndex(0, 0, k));
  return w;
}
}  // namespace

TEST_CASE("b1-sw: Orchestrator reproduces the direct SWE path (L2 < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const double dt = 5.0;
  const long long nsteps = 50;
  const double t_end = dt * static_cast<double>(nsteps);
  const double init_eta = -2.0;

  // --- Direct path (mirrors b1_sw_runner / P4) ---
  const std::vector<double> bath = readDoubles(kB1Dir + "/bath");
  const auto rainflat = readDoubles(kB1Dir + "/rain");
  std::vector<double> rt(rainflat.size() / 2), rv(rainflat.size() / 2);
  for (size_t k = 0; k < rt.size(); ++k) {
    rt[k] = rainflat[2 * k];
    rv[k] = rainflat[2 * k + 1];
  }
  const int gny = static_cast<int>(bath.size());
  double bmin = bath[0];
  for (double b : bath) bmin = std::min(bmin, b);
  const double offset = -bmin;

  MpiComm mc(1, gny, 1, 1);
  Grid grid(mc.localNx(), mc.localNy(), 1, 80.0, 80.0, 0.1);
  SweSolver swe(grid, &mc);
  SweParams p;
  p.gravity = 9.81;
  p.manning = 0.019;
  p.min_depth = 1.0e-8;
  p.visc_x = 1.0e-6;
  p.visc_y = 1.0e-6;
  p.hD = 0.1;
  p.wtfh = 1.0e-8;
  p.dt = dt;
  p.offset = offset;
  swe.setParams(p);
  RealArr1DHost bed("bed", static_cast<size_t>(gny));
  for (int j = 0; j < gny; ++j) bed(static_cast<size_t>(j)) = bath[static_cast<size_t>(j)] + offset;
  swe.setBathymetry(bed);
  swe.initializeState(init_eta);
  Decomp2D dd(mc);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-12;
  PetscLinearSolver solver(cfg);
  swe.attachSolver(solver, dd);
  double tc = 0.0;
  for (long long step = 0; step < nsteps; ++step) {
    tc += dt;
    swe.advanceStep(interpBc(rt, rv, tc), 0.0);
  }
  RealArr1DHost direct = ownedDepth(swe);

  // --- Orchestrator path (same config) ---
  const std::string out = std::string(kTmp) + "/parity_b1.h5";
  Config c = Config::fromString(b1Config(out, init_eta, dt, t_end, nsteps), kB1Dir);
  Orchestrator orch;
  orch.initialize(c);
  orch.run();
  RealArr1DHost driver = ownedDepth(*orch.swe());

  const double l2 = relL2(driver, direct);
  std::fprintf(stderr, "  b1-sw orchestrator-vs-direct rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-10);
}

TEST_CASE("b2-gw: Orchestrator reproduces the direct RE path (L2 < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nz = 50;
  const double dt = 1.0e-4;
  const long long nsteps = 50;
  const double t_end = 1.0e6;  // large; max_steps governs
  const double init_wc = 0.033;

  // --- Direct path (mirrors P5 adaptive run) ---
  MpiComm mc(1, 1, 1, 1);
  Grid grid(1, 1, nz, 1.0, 1.0, 0.01, 1.0);
  ReSolver re(grid, &mc);
  ReParams rp;
  rp.soil.alpha = 1.43;
  rp.soil.n = 1.56;
  rp.soil.theta_s = 0.33;
  rp.soil.theta_r = 0.0;
  rp.soil.Ks_z = 2.89e-6;
  rp.soil.Ss = 1.0e-5;
  rp.soil.use_vg = true;
  rp.dx = 1.0;
  rp.dy = 1.0;
  rp.dz = 0.01;
  rp.botz = -1.0;
  rp.dt = dt;
  rp.dt_min = dt;
  rp.dt_max = 2.0;
  rp.co_max = 2.0;
  rp.adaptive_dt = true;
  rp.use_corrector = true;
  rp.bc_type[5] = 1;
  rp.htop = 0.0;
  re.setParams(rp);
  re.initializeUniformColumn(init_wc);
  Decomp3D dd(mc, nz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "sor";
  cfg.rtol = 1.0e-10;
  PetscLinearSolver solver(cfg);
  re.attachSolver(solver, dd);
  for (long long step = 0; step < nsteps; ++step) re.advanceStep();
  RealArr1DHost direct = wcColumn(re);

  // --- Orchestrator path (same config) ---
  const std::string out = std::string(kTmp) + "/parity_b2.h5";
  Config c = Config::fromString(gwConfig(out, nz, dt, t_end, nsteps, init_wc), "");
  Orchestrator orch;
  orch.initialize(c);
  orch.run();
  RealArr1DHost driver = wcColumn(*orch.re());

  const double l2 = relL2(driver, direct);
  std::fprintf(stderr, "  b2-gw orchestrator-vs-direct rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-10);
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
