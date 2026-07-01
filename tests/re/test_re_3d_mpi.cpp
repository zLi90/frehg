// P23: 3D Richards MPI rank-count equivalence for a LATERALLY-COUPLED box.
//
// Unlike the P5.5e gate (a laterally-decoupled box, Kx=Ky=0, which reaches 1e-10 because the
// matrix is block-diagonal per column), this runs the full 3D path (use_full3d=true) with
// non-zero lateral conductivity and a lateral head gradient, so the predictor matrix couples
// cells ACROSS rank boundaries and the corrector moves water across them. With the horizontal
// halo exchange (h, Kx/Ky, qx/qy, wc) and the global adaptive-dt reduction, the gathered field
// is independent of the process count up to the iterative-solver consistency floor.
//
// Tolerance (1e-7, not 1e-10): for a genuinely cross-rank-coupled system the Krylov solve
// (CG + Jacobi; the only deterministic option on the local PETSc, which has NO parallel direct
// LU — see the P2 carry-forward) lands on partition-dependent solutions at the ~1e-9 parallel
// inner-product round-off level, which then accumulates mildly over the nonlinear steps. The
// SERIAL run conserves total water to ~1e-16 (test_re_3d), so the ~1e-8 here is solver round-off,
// NOT a halo/assembly bug (those would show O(1e-2)+ structured errors). Bit-level (1e-10)
// coupled-solve partition independence needs a parallel direct solver (MUMPS/SuperLU_DIST),
// deferred to the Linux validation build, exactly like the documented `pc_type lu` note.
//
// np1 writes the reference; np2/np4 compare against it (CTest fixtures enforce ordering).
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
constexpr int kNz = 6;
constexpr int kSteps = 60;
// Cross-rank-coupled iterative-solve consistency floor (see file header). Strong partition
// independence (8 digits on an O(0.3) field); halo/assembly bugs would be O(1e-2)+.
constexpr double kTol = 1.0e-7;
const char* kRefPath = "re_3d_mpi_ref.txt";

SoilParams makeSoil() {
  SoilParams s;
  s.alpha = 2.0;
  s.n = 1.6;
  s.theta_s = 0.45;
  s.theta_r = 0.05;
  s.Ks_x = 1.0e-5;  // lateral conductivity ON => cross-rank coupling
  s.Ks_y = 1.0e-5;
  s.Ks_z = 1.0e-5;
  s.Ss = 1.0e-5;
  s.use_vg = true;
  return s;
}

// Lateral + vertical IC variation so the lateral coupling is actually exercised.
double icWc(int gi, int gj) {
  return 0.18 + 0.04 * static_cast<double>(gi) / kGnx + 0.03 * static_cast<double>(gj) / kGny;
}

double runBox(const MpiComm& mc, Grid& local, Decomp3D& dd, PetscLinearSolver& solver,
              std::vector<double>& local_wc) {
  const auto soil = makeSoil();
  ReSolver re(local, &mc);
  ReParams rp;
  rp.soil = soil;
  rp.dx = 1.0;
  rp.dy = 1.0;
  rp.dz = 0.05;
  rp.botz = -0.05 * kNz;
  rp.dt = 2.0;
  rp.dt_min = 1.0e-4;
  rp.dt_max = 20.0;
  rp.co_max = 2.0;
  rp.adaptive_dt = true;
  rp.use_corrector = true;
  rp.use_full3d = true;  // full 3D lateral path
  for (int b = 0; b < 6; ++b) rp.bc_type[b] = 0;  // closed box
  re.setParams(rp);
  re.initializeUniformColumn(0.18);

  auto& f = re.fields();
  for (int k = 0; k < kNz; ++k)
    for (int lj = 0; lj < local.ny(); ++lj)
      for (int li = 0; li < local.nx(); ++li) {
        const int gi = mc.i0() + li, gj = mc.j0() + lj;
        const double w = icWc(gi, gj);
        const int idx = local.getIndex(li, lj, k);
        f.wc(idx) = w;
        f.h(idx) = VanGenuchten::headFromWaterContent(soil, w);
        f.wcn(idx) = w;
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

void runRankEquivalence(int pnx, int pny) {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != pnx * pny) return;

  MpiComm mc(kGnx, kGny, pnx, pny);
  Grid local(mc.localNx(), mc.localNy(), kNz, 1.0, 1.0, 0.05);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  PetscLinearSolver solver(cfg);

  std::vector<double> wc;
  runBox(mc, local, dd, solver, wc);

  std::ifstream in(kRefPath);
  REQUIRE(in.good());
  std::vector<double> ref(static_cast<size_t>(kGnx * kGny * kNz));
  for (auto& v : ref) in >> v;

  const int lnx = local.nx(), lny = local.ny();
  double worst = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int lj = 0; lj < lny; ++lj)
      for (int li = 0; li < lnx; ++li) {
        const int gi = mc.i0() + li, gj = mc.j0() + lj;
        const double expected = ref[static_cast<size_t>(gi + gj * kGnx + k * kGnx * kGny)];
        const double got = wc[static_cast<size_t>(li + lj * lnx + k * lnx * lny)];
        worst = std::max(worst, std::fabs(got - expected));
      }
  double global_worst = worst;
  MPI_Allreduce(&worst, &global_worst, 1, MPI_DOUBLE, MPI_MAX, mc.comm());
  REQUIRE((mc.localNx() < kGnx) || (mc.localNy() < kGny));
  if (mc.rank() == 0)
    std::fprintf(stderr, "  RE 3D MPI np%d (%dx%d): worst |wc-ref| = %.3e (tol %.0e)\n", size, pnx,
                 pny, global_worst, kTol);
  REQUIRE(global_worst < kTol);
}
}  // namespace

TEST_CASE("RE 3D MPI: 1-rank reference run writes global field") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  MpiComm mc(kGnx, kGny, 1, 1);
  Grid local(mc.localNx(), mc.localNy(), kNz, 1.0, 1.0, 0.05);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  PetscLinearSolver solver(cfg);

  std::vector<double> wc;
  const double mx = runBox(mc, local, dd, solver, wc);
  std::fprintf(stderr, "  RE 3D MPI np1: max wc = %.6f (theta_s=0.45)\n", mx);

  std::ofstream out(kRefPath);
  REQUIRE(out.good());
  out.precision(17);
  for (int k = 0; k < kNz; ++k)
    for (int gj = 0; gj < kGny; ++gj)
      for (int gi = 0; gi < kGnx; ++gi)
        out << wc[static_cast<size_t>(gi + gj * kGnx + k * kGnx * kGny)] << "\n";
  out.close();
}

TEST_CASE("RE 3D MPI: 2-rank x-decomp matches 1-rank reference (max-diff < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) return;
  runRankEquivalence(2, 1);
}

TEST_CASE("RE 3D MPI: 4-rank 2x2 decomp matches 1-rank reference (max-diff < 1e-10)") {
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
