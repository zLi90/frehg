// P5.5e: RE solver MPI rank-count equivalence gate.
//
// b2-gw is intrinsically a 1x1 column (no horizontal cells to decompose; z is on-rank),
// so it cannot exercise horizontal decomposition. The implemented RE corrector is
// vertical (qx=qy=0); this gate therefore runs a multi-COLUMN box with no-flow lateral
// boundaries (every column is vertically coupled but laterally independent: Kx=Ky=0 for
// unsaturated cells), distributed across ranks. With the global adaptive-dt reduction the
// dt sequence is identical on every rank, so the gathered field must be independent of the
// process count. Gate: max abs difference vs the 1-rank reference < 1e-10.
//
// The 1-rank run (test_re_mpi, np1) writes the reference; the np2/np4 runs compare against
// it (CTest DEPENDS enforces ordering).
#include <petscksp.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "re/ReSolver.hpp"
#include "re/VanGenuchten.hpp"

using namespace frehg2;

namespace {
constexpr int kGnx = 4;
constexpr int kGny = 4;
constexpr int kNz = 30;
constexpr int kSteps = 300;
const char* kRefPath = "re_mpi_ref.txt";

SoilParams makeSoil() {
  SoilParams soil;
  soil.alpha = 1.43;
  soil.n = 1.56;
  soil.theta_s = 0.33;
  soil.theta_r = 0.0;
  soil.Ks_z = 2.89e-6;
  soil.Ks_x = 2.89e-6;
  soil.Ks_y = 2.89e-6;
  soil.Ss = 1.0e-5;
  soil.use_vg = true;
  return soil;
}

// Spatially varied, always-unsaturated IC so lateral conductivity stays zero.
double icWaterContent(int gi, int gj) {
  return 0.08 + 0.01 * static_cast<double>(gi) + 0.005 * static_cast<double>(gj);
}

// Run the box on this rank's subdomain for kSteps; fill local_wc[li + lj*lnx + k*lnx*lny]
// and report the local max wc.
double runBox(const MpiComm& mc, Grid& local, Decomp3D& dd, PetscLinearSolver& solver,
              std::vector<double>& local_wc) {
  const auto soil = makeSoil();
  ReSolver re(local, &mc);
  ReParams rp;
  rp.soil = soil;
  rp.dx = 1.0;
  rp.dy = 1.0;
  rp.dz = 0.01;
  rp.botz = -0.30;
  rp.dt = 1.0e-4;
  rp.dt_min = 1.0e-4;
  rp.dt_max = 2.0;
  rp.co_max = 2.0;
  rp.adaptive_dt = true;
  rp.use_corrector = true;
  rp.use_full3d = false;
  for (int b = 0; b < 6; ++b) rp.bc_type[b] = 0;  // all no-flow
  re.setParams(rp);
  re.initializeUniformColumn(0.08);

  // Override with the varied IC by global column index.
  auto& f = re.fields();
  for (int k = 0; k < kNz; ++k)
    for (int lj = 0; lj < local.ny(); ++lj)
      for (int li = 0; li < local.nx(); ++li) {
        const int gi = mc.i0() + li;
        const int gj = mc.j0() + lj;
        const double wc0 = icWaterContent(gi, gj);
        const int idx = local.getIndex(li, lj, k);
        f.wc(idx) = wc0;
        f.h(idx) = VanGenuchten::headFromWaterContent(soil, wc0);
        f.wcn(idx) = wc0;
        f.hn(idx) = f.h(idx);
      }

  re.attachSolver(solver, dd);
  for (int step = 0; step < kSteps; ++step) re.advanceStep();

  const int lnx = local.nx(), lny = local.ny();
  local_wc.assign(static_cast<size_t>(lnx * lny * kNz), 0.0);
  double m = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int lj = 0; lj < lny; ++lj)
      for (int li = 0; li < lnx; ++li) {
        const double v = f.wc(local.getIndex(li, lj, k));
        local_wc[static_cast<size_t>(li + lj * lnx + k * lnx * lny)] = v;
        m = std::max(m, v);
      }
  return m;
}
}  // namespace

TEST_CASE("RE MPI: 1-rank reference run writes global field") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  MpiComm mc(kGnx, kGny, 1, 1);
  Grid local(mc.localNx(), mc.localNy(), kNz, 1.0, 1.0, 0.01);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  PetscLinearSolver solver(cfg);

  std::vector<double> wc;
  const double mx = runBox(mc, local, dd, solver, wc);  // size 1: local == global

  // Guard the laterally-independent assumption: nothing saturated.
  std::fprintf(stderr, "  RE MPI np1: max wc = %.6f (theta_s=0.33)\n", mx);
  REQUIRE(mx < 0.33);

  std::ofstream out(kRefPath);
  REQUIRE(out.good());
  out.precision(17);
  for (int k = 0; k < kNz; ++k)
    for (int gj = 0; gj < kGny; ++gj)
      for (int gi = 0; gi < kGnx; ++gi)
        out << wc[static_cast<size_t>(gi + gj * kGnx + k * kGnx * kGny)] << "\n";
  out.close();
}

void runRankEquivalence(int pnx, int pny);
void runRankEquivalence(int pnx, int pny) {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != pnx * pny) return;

  MpiComm mc(kGnx, kGny, pnx, pny);
  Grid local(mc.localNx(), mc.localNy(), kNz, 1.0, 1.0, 0.01);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  PetscLinearSolver solver(cfg);

  std::vector<double> wc;
  runBox(mc, local, dd, solver, wc);

  // Load the global reference field.
  std::ifstream in(kRefPath);
  REQUIRE(in.good());
  std::vector<double> ref(static_cast<size_t>(kGnx * kGny * kNz));
  for (auto& v : ref) in >> v;

  const int lnx = local.nx(), lny = local.ny();
  double worst = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int lj = 0; lj < lny; ++lj)
      for (int li = 0; li < lnx; ++li) {
        const int gi = mc.i0() + li;
        const int gj = mc.j0() + lj;
        const double expected = ref[static_cast<size_t>(gi + gj * kGnx + k * kGnx * kGny)];
        const double got = wc[static_cast<size_t>(li + lj * lnx + k * lnx * lny)];
        worst = std::max(worst, std::fabs(got - expected));
      }
  double global_worst = worst;
  MPI_Allreduce(&worst, &global_worst, 1, MPI_DOUBLE, MPI_MAX, mc.comm());
  // Confirm the domain actually decomposed on this rank layout.
  REQUIRE((mc.localNx() < kGnx) || (mc.localNy() < kGny));
  if (mc.rank() == 0)
    std::fprintf(stderr, "  RE MPI np%d (%dx%d): worst |wc-ref| = %.3e\n", size, pnx, pny,
                 global_worst);
  REQUIRE(global_worst < 1.0e-10);
}

TEST_CASE("RE MPI: 2-rank x-decomp matches 1-rank reference (max-diff < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) return;
  runRankEquivalence(2, 1);
}

TEST_CASE("RE MPI: 4-rank 2x2 decomp matches 1-rank reference (max-diff < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 4) return;
  runRankEquivalence(2, 2);
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
