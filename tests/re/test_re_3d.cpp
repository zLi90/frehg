// P23: 3D Richards solver (lateral flux) + fully heterogeneous per-cell soil — serial gates.
//
//   1. Closed-box conservation: a closed (all no-flux) box with Ss=0 and an unsaturated
//      heterogeneous IC conserves total water (Sum wc*V) to machine precision while the
//      predictor-corrector redistributes it in 3D (dx != dy exercises the face-area factor).
//   2. Lateral activity: turning lateral conductivity on (use_full3d) produces a DIFFERENT
//      final field than the vertical-only path, proving qx/qy actually move water.
//   3. Cell-by-cell soil: a single-class 3D SoilMap reproduces the uniform run bit-for-bit;
//      a vertically layered 2-class map produces a different (per-cell) result.
#include <petscksp.h>

#include <cmath>
#include <cstdio>
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
#include "soil/SoilMap.hpp"

using namespace frehg2;

namespace {
constexpr int kNx = 4;
constexpr int kNy = 4;
constexpr int kNz = 8;
constexpr int kSteps = 40;
constexpr double kDx = 2.0;  // dx != dy on purpose (exercises Ax=dy*dz vs Ay=dx*dz)
constexpr double kDy = 1.0;
constexpr double kDz = 0.05;

SoilParams makeSoil(double Ksx, double Ksy, double Ksz) {
  SoilParams s;
  s.alpha = 2.0;
  s.n = 1.6;
  s.theta_s = 0.45;
  s.theta_r = 0.05;
  s.Ks_x = Ksx;
  s.Ks_y = Ksy;
  s.Ks_z = Ksz;
  s.Ss = 0.0;  // exact volume conservation in the corrector
  s.use_vg = true;
  return s;
}

// Smoothly varying, always-unsaturated IC so neither clamp nor the realloc saturation branch
// fires (the flux divergence then conserves Sum wc*V exactly).
double icWc(int gi, int gj, int /*k*/) {
  return 0.20 + 0.03 * std::sin(0.9 * gi) + 0.02 * std::cos(0.7 * gj);
}

ReParams baseParams(const SoilParams& soil, bool full3d) {
  ReParams rp;
  rp.soil = soil;
  rp.dx = kDx;
  rp.dy = kDy;
  rp.dz = kDz;
  rp.botz = -kDz * kNz;
  rp.dt = 5.0;
  rp.adaptive_dt = false;  // fixed dt: a clean, reproducible conservation check
  rp.use_corrector = true;
  rp.use_full3d = full3d;
  for (int b = 0; b < 6; ++b) rp.bc_type[b] = 0;  // closed box (no flux on all six faces)
  return rp;
}

void setIc(ReSolver& re, const MpiComm& mc, const SoilParams& soil) {
  auto& f = re.fields();
  for (int k = 0; k < kNz; ++k)
    for (int lj = 0; lj < mc.localNy(); ++lj)
      for (int li = 0; li < mc.localNx(); ++li) {
        const int gi = mc.i0() + li, gj = mc.j0() + lj;
        const double w = icWc(gi, gj, k);
        const int idx = re.grid().getIndex(li, lj, k);
        f.wc(idx) = w;
        f.h(idx) = VanGenuchten::headFromWaterContent(soil, w);
        f.wcn(idx) = w;
        f.hn(idx) = f.h(idx);
      }
}

double totalWater(ReSolver& re) {
  auto& f = re.fields();
  const double V = kDx * kDy * kDz;
  double sum = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int j = 0; j < kNy; ++j)
      for (int i = 0; i < kNx; ++i) sum += f.wc(re.grid().getIndex(i, j, k)) * V;
  return sum;
}

// Run a serial 3D box; optionally attach a SoilMap. Returns the final wc field (i+j*nx+k*nx*ny).
std::vector<double> run3d(bool full3d, const SoilParams& soil, const SoilMap* map,
                          double& total_initial, double& total_final) {
  MpiComm mc(kNx, kNy, 1, 1);
  Grid local(kNx, kNy, kNz, kDx, kDy, kDz);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  PetscLinearSolver solver(cfg);

  ReSolver re(local, &mc);
  re.setParams(baseParams(soil, full3d));
  if (map != nullptr) re.setSoilMap(map);
  re.initializeGeometry();
  setIc(re, mc, soil);
  re.finalizeInitialState();
  re.attachSolver(solver, dd);

  total_initial = totalWater(re);
  for (int s = 0; s < kSteps; ++s) re.advanceStep();
  total_final = totalWater(re);

  std::vector<double> out(static_cast<size_t>(kNx * kNy * kNz));
  auto& f = re.fields();
  for (int k = 0; k < kNz; ++k)
    for (int j = 0; j < kNy; ++j)
      for (int i = 0; i < kNx; ++i)
        out[static_cast<size_t>(i + j * kNx + k * kNx * kNy)] = f.wc(local.getIndex(i, j, k));
  return out;
}
}  // namespace

TEST_CASE("RE 3D: closed heterogeneous box conserves total water to machine precision") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const SoilParams soil = makeSoil(1.0e-5, 1.0e-5, 1.0e-5);
  double t0 = 0.0, t1 = 0.0;
  const std::vector<double> wc3d = run3d(true, soil, nullptr, t0, t1);

  const double rel = std::fabs(t1 - t0) / t0;
  std::fprintf(stderr, "  RE 3D conservation: initial=%.12e final=%.12e rel=%.3e\n", t0, t1, rel);
  REQUIRE(rel < 1.0e-10);
  // State must remain physical (bounded by [theta_r, theta_s]).
  for (double w : wc3d) {
    REQUIRE(w >= soil.theta_r - 1e-12);
    REQUIRE(w <= soil.theta_s + 1e-12);
  }
}

TEST_CASE("RE 3D: lateral conductivity changes the result vs the vertical-only path") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  double a0, a1, b0, b1;
  const std::vector<double> with_lat = run3d(true, makeSoil(1.0e-5, 1.0e-5, 1.0e-5), nullptr, a0, a1);
  // Lateral conductivity zero => qx=qy=0 => purely vertical (the P5 behaviour).
  const std::vector<double> no_lat = run3d(true, makeSoil(0.0, 0.0, 1.0e-5), nullptr, b0, b1);

  double worst = 0.0;
  for (size_t i = 0; i < with_lat.size(); ++i)
    worst = std::max(worst, std::fabs(with_lat[i] - no_lat[i]));
  std::fprintf(stderr, "  RE 3D lateral activity: max|wc_3d - wc_vert| = %.3e\n", worst);
  REQUIRE(worst > 1.0e-7);   // lateral flux measurably redistributes water (>> round-off)
  REQUIRE(std::fabs(a1 - a0) / a0 < 1.0e-10);  // both paths still conserve
  REQUIRE(std::fabs(b1 - b0) / b0 < 1.0e-10);
}

TEST_CASE("RE 3D: single-class 3D soil map reproduces the uniform run bit-for-bit") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const SoilParams soil = makeSoil(1.0e-5, 1.0e-5, 1.0e-5);
  double u0, u1, m0, m1;
  const std::vector<double> uniform = run3d(true, soil, nullptr, u0, u1);

  SoilMap one;
  one.setClasses({soil});
  one.setClassIndex3D(kNx, kNy, kNz, std::vector<int>(static_cast<size_t>(kNx * kNy * kNz), 0));
  const std::vector<double> mapped = run3d(true, soil, &one, m0, m1);

  double worst = 0.0;
  for (size_t i = 0; i < uniform.size(); ++i)
    worst = std::max(worst, std::fabs(uniform[i] - mapped[i]));
  std::fprintf(stderr, "  RE 3D single-class map vs uniform: max|diff| = %.3e\n", worst);
  REQUIRE(worst == 0.0);  // bit-identical
}

TEST_CASE("RE 3D: per-cell (vertically layered) soil map changes the result") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // Two classes: fast on top half, slow on the bottom half (a per-cell vertical contrast).
  const SoilParams fast = makeSoil(5.0e-5, 5.0e-5, 5.0e-5);
  const SoilParams slow = makeSoil(2.0e-6, 2.0e-6, 2.0e-6);
  SoilMap layered;
  layered.setClasses({fast, slow});
  std::vector<int> idx(static_cast<size_t>(kNx * kNy * kNz), 0);
  for (int k = 0; k < kNz; ++k)
    for (int j = 0; j < kNy; ++j)
      for (int i = 0; i < kNx; ++i)
        idx[static_cast<size_t>(i + j * kNx + k * kNx * kNy)] = (k < kNz / 2) ? 0 : 1;
  layered.setClassIndex3D(kNx, kNy, kNz, std::move(idx));

  double a0, a1, b0, b1;
  const std::vector<double> all_fast = run3d(true, fast, nullptr, a0, a1);
  const std::vector<double> mixed = run3d(true, fast, &layered, b0, b1);

  double worst = 0.0;
  for (size_t i = 0; i < all_fast.size(); ++i)
    worst = std::max(worst, std::fabs(all_fast[i] - mixed[i]));
  std::fprintf(stderr, "  RE 3D layered vs all-fast: max|diff| = %.3e\n", worst);
  REQUIRE(worst > 1.0e-6);  // the per-cell soil contrast changes the solution
}

TEST_CASE("RE bottom fixed-flux (qbot, bc_type[4]==2) adds exactly the prescribed mass") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const SoilParams soil = makeSoil(1.0e-5, 1.0e-5, 1.0e-5);
  const double qbot = 1.0e-7;  // m/s, >0 = upward inflow INTO the domain at the bottom
  const double dt = 5.0;
  const int steps = 40;

  MpiComm mc(kNx, kNy, 1, 1);
  Grid local(kNx, kNy, kNz, kDx, kDy, kDz);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-13;
  PetscLinearSolver solver(cfg);

  ReParams rp = baseParams(soil, /*full3d=*/true);
  rp.dt = dt;
  rp.bc_type[4] = 2;  // bottom (z+) fixed flux
  rp.qbot = qbot;

  ReSolver re(local, &mc);
  re.setParams(rp);
  re.initializeGeometry();
  setIc(re, mc, soil);
  re.finalizeInitialState();
  re.attachSolver(solver, dd);

  const double t0 = totalWater(re);
  for (int s = 0; s < steps; ++s) re.advanceStep();
  const double t1 = totalWater(re);

  // Prescribed inflow volume = qbot * (dx*dy) per column, per step, over all columns/steps.
  const double expected = qbot * kDx * kDy * dt * static_cast<double>(steps) *
                          static_cast<double>(kNx * kNy);
  const double added = t1 - t0;
  const double rel = std::fabs(added - expected) / expected;
  std::fprintf(stderr, "  RE qbot mass balance: added=%.6e expected=%.6e rel=%.3e\n", added,
               expected, rel);
  REQUIRE(added > 0.0);          // inflow actually adds water (not a silent no-op)
  REQUIRE(rel < 1.0e-6);         // conserved to the prescribed flux
  // State stays physical and bounded.
  auto& f = re.fields();
  for (int k = 0; k < kNz; ++k)
    for (int j = 0; j < kNy; ++j)
      for (int i = 0; i < kNx; ++i) {
        const double w = f.wc(local.getIndex(i, j, k));
        REQUIRE(std::isfinite(w));
        REQUIRE(w >= soil.theta_r - 1e-12);
        REQUIRE(w <= soil.theta_s + 1e-12);
      }
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rc = frehg2test::runAll();
  PetscFinalize();
  MPI_Finalize();
  Kokkos::finalize();
  return rc;
}
