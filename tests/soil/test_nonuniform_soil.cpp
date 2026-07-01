// P13 acceptance gate (manufactured, user-authorized 2026-06-28): non-uniform soil correctness.
//
// The plan's named gate (b4-govindaraju "L2<1e-3 vs legacy multi-class soil") cannot be used:
// the registered b4-govindaraju is a SERGHEI overland-flow (SWE) benchmark with NO soil/GW
// (DEM + rainfall + Chezy friction, InfModel:none). Per .cursorrules the reference_registry is
// authoritative and tolerances must not be relabeled, so P13 is gated by this self-consistent
// numerical test instead (plus the b2-gw uniform bit-identical regression).
//
// Construction: a multi-column box with two soil classes (left half = class A, right half =
// class B, differing only in saturated conductivity). All lateral boundaries are no-flow and
// every cell stays unsaturated, so computeKFace() zeroes Kx/Ky and the columns are laterally
// INDEPENDENT (the global matrix is block-diagonal per column). With a fixed dt, gravity
// redistributes water within each column at a rate set by that column's class.
//
// Two rigorous assertions:
//   (1) Per-cell dispatch: every class-A column of the NON-UNIFORM run matches the same column
//       of a UNIFORM class-A run, and every class-B column matches a UNIFORM class-B run, to
//       < 1e-9 (the only difference is PETSc CG cross-block coupling at rtol=1e-14). This proves
//       the SoilMap is read per cell AND that the uniform single-class path is unchanged (P5).
//   (2) Meaningfulness: the uniform A and B runs differ by a wide margin (> 1e-3) at matching
//       columns, so the per-cell soil genuinely changes the solution (a no-op map would fail).
#include <petscksp.h>

#include <algorithm>
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
#include "soil/SoilMap.hpp"

using namespace frehg2;

namespace {
constexpr int kNx = 6;   // columns 0..2 = class A, 3..5 = class B
constexpr int kNy = 1;
constexpr int kNz = 24;
constexpr int kSteps = 150;
constexpr double kDx = 1.0, kDy = 1.0, kDz = 0.01;
constexpr double kIcWc = 0.22;  // < theta_s (keeps cells unsaturated => Kx=Ky=0)

SoilParams classWith(double ks) {
  SoilParams s;
  s.alpha = 1.43;
  s.n = 1.56;
  s.theta_s = 0.33;
  s.theta_r = 0.0;
  s.Ks_x = ks;
  s.Ks_y = ks;
  s.Ks_z = ks;
  s.Ss = 1.0e-5;
  s.use_vg = true;
  return s;
}
const SoilParams kClassA = classWith(1.0e-5);
const SoilParams kClassB = classWith(1.0e-4);  // 10x faster redistribution

ReParams baseParams(const SoilParams& uniform_class) {
  ReParams rp;
  rp.soil = uniform_class;
  rp.dx = kDx;
  rp.dy = kDy;
  rp.dz = kDz;
  rp.botz = -kDz * kNz;
  rp.dt = 2.0;          // fixed dt => identical time sequence across all three runs
  rp.dt_min = 2.0;
  rp.dt_max = 2.0;
  rp.adaptive_dt = false;
  rp.use_corrector = true;
  rp.use_full3d = false;
  for (int b = 0; b < 6; ++b) rp.bc_type[b] = 0;  // all no-flow
  return rp;
}

// Run a kNx x kNy x kNz box for kSteps. If `map` is non-null it is attached (non-uniform soil),
// otherwise the uniform rp.soil path is used. Returns the full wc field (col i + k*kNx) and the
// max wc seen (to guard the unsaturated/lateral-decoupling assumption).
double runBox(const SoilParams& uniform_class, const SoilMap* map, std::vector<double>& wc) {
  MpiComm mc(kNx, kNy, 1, 1);
  Grid grid(mc.localNx(), mc.localNy(), kNz, kDx, kDy, kDz);
  Decomp3D dd(mc, kNz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-14;
  PetscLinearSolver solver(cfg);

  ReSolver re(grid, &mc);
  re.setParams(baseParams(uniform_class));
  if (map != nullptr) re.setSoilMap(map);
  re.initializeUniformColumn(kIcWc);
  re.attachSolver(solver, dd);
  for (int step = 0; step < kSteps; ++step) re.advanceStep();

  const auto& f = re.fields();
  wc.assign(static_cast<size_t>(kNx) * kNz, 0.0);
  double mx = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int i = 0; i < kNx; ++i) {
      const double v = f.wc(grid.getIndex(i, 0, k));
      wc[static_cast<size_t>(i + k * kNx)] = v;
      mx = std::max(mx, v);
    }
  return mx;
}
}  // namespace

TEST_CASE("non-uniform soil: per-cell SoilMap dispatch matches per-class uniform runs") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // single-rank gate

  // Non-uniform run: left half class A, right half class B.
  SoilMap map;
  map.setClasses({kClassA, kClassB});
  std::vector<int> idx(static_cast<size_t>(kNx) * kNy, 0);
  for (int i = 0; i < kNx; ++i) idx[static_cast<size_t>(i)] = (i < kNx / 2) ? 0 : 1;
  map.setClassIndex(kNx, kNy, idx);

  std::vector<double> wc_nu, wc_a, wc_b;
  const double mx_nu = runBox(kClassA, &map, wc_nu);
  const double mx_a = runBox(kClassA, nullptr, wc_a);
  const double mx_b = runBox(kClassB, nullptr, wc_b);

  // Lateral-decoupling guard: every cell stayed unsaturated (theta_s = 0.33).
  std::fprintf(stderr, "  non-uniform soil: max wc  nu=%.4f a=%.4f b=%.4f (theta_s=0.33)\n",
               mx_nu, mx_a, mx_b);
  REQUIRE(mx_nu < 0.33);
  REQUIRE(mx_a < 0.33);
  REQUIRE(mx_b < 0.33);

  // (1) Per-cell dispatch: NU class-A columns == uniform-A; NU class-B columns == uniform-B.
  double worst_match = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int i = 0; i < kNx; ++i) {
      const double got = wc_nu[static_cast<size_t>(i + k * kNx)];
      const double ref =
          (i < kNx / 2) ? wc_a[static_cast<size_t>(i + k * kNx)]
                        : wc_b[static_cast<size_t>(i + k * kNx)];
      worst_match = std::max(worst_match, std::fabs(got - ref));
    }
  std::fprintf(stderr, "  non-uniform soil: worst |nu - per-class uniform| = %.3e\n", worst_match);
  REQUIRE(worst_match < 1.0e-9);

  // (2) Meaningfulness: class A and class B genuinely differ at matching columns (a no-op or
  // ignored soil map would make A == B and silently pass assertion (1)).
  double margin = 0.0;
  for (int k = 0; k < kNz; ++k)
    for (int i = 0; i < kNx; ++i)
      margin = std::max(margin, std::fabs(wc_a[static_cast<size_t>(i + k * kNx)] -
                                          wc_b[static_cast<size_t>(i + k * kNx)]));
  std::fprintf(stderr, "  non-uniform soil: max |uniform_A - uniform_B| = %.3e\n", margin);
  REQUIRE(margin > 1.0e-3);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  const int rc = frehg2test::runAll();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  MPI_Finalize();
  Kokkos::finalize();
  return grc;
}
