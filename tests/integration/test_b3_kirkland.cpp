// b3-kirkland FULL 2-D layered-soil infiltration (Kirkland et al. 1992), SERGHEI suite.
//
// This is the full 2-D (x-z) port enabled by P23 (3-D Richards lateral flux + per-cell soil) plus
// the new fixed-flux top BC (legacy bctype_GW[5]==2 / qtop, wired as groundwater.recharge). The
// solver is the sanctioned PCA predictor-corrector (NO Picard/Newton). The reference is the
// digitized Kirkland h=0 / h=-400 hydraulic-head contours embedded in the SERGHEI makeplot.py
// (ref0 / ref400). The registry gates b3 `review` (tolerance: null), so this asserts review-tier
// physics: stable + bounded state, mass conservation of the prescribed recharge, a saturated plume
// under the source, and a contour match to the Kirkland reference within a loose band (the PCA
// scheme + arithmetic K-mean give an approximate, not bit-exact, match — Picard would be needed
// for tighter agreement, and it is excluded by policy).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <mpi.h>
#include <string>
#include <vector>

#include "core/Grid.hpp"
#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/benchmark_util.hpp"
#include "re/ReSolver.hpp"
#include "soil/SoilMap.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;
const char* kSrc = FREHG2_SOURCE_DIR;

constexpr int NX = 50, NY = 1, NZ = 30;
constexpr double DX = 0.1, HGT = NZ * DX;

// Digitized Kirkland reference contours from legacy/benchmarks/b3-kirkland/makeplot.py:
// (x [m], z-above-bottom [m]) pairs for the h=0 and h=-400 hydraulic-head contours.
const std::array<std::array<double, 2>, 13> kRef0 = {{{1.0881, 2.1851},
                                                      {1.5132, 2.1943},
                                                      {2.0205, 2.2198},
                                                      {2.6435, 2.2246},
                                                      {3.2119, 2.1106},
                                                      {3.7317, 2.1741},
                                                      {3.9631, 2.0063},
                                                      {3.5125, 1.7818},
                                                      {2.9268, 1.7478},
                                                      {2.3033, 1.7440},
                                                      {1.7238, 1.7632},
                                                      {1.2827, 1.8152},
                                                      {1.0151, 1.9801}}};
const std::array<std::array<double, 2>, 13> kRef400 = {{{0.4416, 2.9334},
                                                        {0.4739, 2.4642},
                                                        {0.5570, 1.9486},
                                                        {0.7014, 1.7254},
                                                        {1.0389, 1.4954},
                                                        {1.5459, 1.3936},
                                                        {2.0466, 1.3551},
                                                        {3.0030, 1.3560},
                                                        {3.8003, 1.4410},
                                                        {4.2949, 1.7139},
                                                        {4.4776, 1.9707},
                                                        {4.5197, 2.4166},
                                                        {4.5554, 2.8034}}};

// Deepest z-above-bottom at column i where head >= level (the wetting/contour front), or NaN.
double frontZ(const std::vector<double>& hxz, int i, double level) {
  double zf = std::numeric_limits<double>::quiet_NaN();
  for (int k = 0; k < NZ; ++k) {
    if (hxz[static_cast<size_t>(k * NX + i)] >= level) zf = HGT - (k + 0.5) * DX;
  }
  return zf;
}

// RMS |Δz| of the simulated `level` front vs the reference points (skips columns with no front).
template <size_t N>
double contourRms(const std::vector<double>& hxz, const std::array<std::array<double, 2>, N>& ref,
                  double level, int& n_used) {
  double sse = 0.0;
  n_used = 0;
  for (const auto& p : ref) {
    int i = static_cast<int>(std::lround(p[0] / DX - 0.5));
    i = std::max(0, std::min(NX - 1, i));
    const double zf = frontZ(hxz, i, level);
    if (std::isfinite(zf)) {
      sse += (zf - p[1]) * (zf - p[1]);
      ++n_used;
    }
  }
  return n_used > 0 ? std::sqrt(sse / n_used) : std::nan("");
}

double gwStorage(const ReSolver& re) {
  const Grid& g = re.grid();
  const auto& f = re.fields();
  const double cellV = g.dx() * g.dy() * g.dz();  // uniform dz here
  double v = 0.0;
  for (int k = 0; k < g.nz(); ++k)
    for (int j = 0; j < g.ny(); ++j)
      for (int i = 0; i < g.nx(); ++i) v += f.wc(g.getIndex(i, j, k)) * cellV;
  return v;
}
}  // namespace

TEST_CASE("b3-kirkland full 2-D layered infiltration: stable, conservative, matches Kirkland") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // serial reference run (the 2-D PCA solve is validated on one rank)

  const std::string dir = std::string(kSrc) + "/benchmarks/b3-kirkland";
  const std::string yaml = dir + "/b3-kirkland.yaml";
  // Full run to t_end = 86400 s (the reference contours are at the final time); ~8.7k adaptive
  // steps, ~1 min. max_steps = 0 => uncapped.
  Config cfg = bench::loadCapped(yaml, dir, 0, std::string(kTmp) + "/b3.h5");

  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.re() != nullptr);
  const ReSolver& re = *orch.re();
  REQUIRE(re.grid().nx() == NX);
  REQUIRE(re.grid().nz() == NZ);
  REQUIRE(re.params().use_full3d);                 // 3-D lateral corrector active
  REQUIRE(re.params().bc_type[5] == 2);            // fixed-flux top BC active
  REQUIRE(re.soilMap() != nullptr);                // per-cell 3-D soil map attached
  REQUIRE(re.soilMap()->is3D());

  const double store0 = gwStorage(re);
  orch.run();

  const Grid& g = re.grid();
  const auto& gf = re.fields();
  const SoilMap& sm = *re.soilMap();

  // (1) Stable + bounded: finite head, wc within each cell's own [theta_r, theta_s].
  std::vector<double> hxz(static_cast<size_t>(NX) * NZ, 0.0);
  for (int k = 0; k < NZ; ++k)
    for (int i = 0; i < NX; ++i) {
      const int c = g.getIndex(i, 0, k);
      const double h = gf.h(c), w = gf.wc(c);
      const SoilParams& sp = sm.classAt(i, 0, k);
      REQUIRE(std::isfinite(h));
      REQUIRE(std::isfinite(w));
      REQUIRE(w >= sp.theta_r - 1.0e-9);
      REQUIRE(w <= sp.theta_s + 1.0e-9);
      hxz[static_cast<size_t>(k * NX + i)] = h;
    }

  // (2) Mass conservation: the storage increase equals the prescribed top recharge volume
  // (|rate| * recharge_area * t_end) to within a few percent (the simplified saturation
  // reallocation and residual clamp account for the small remainder).
  const double store1 = gwStorage(re);
  const double rate = 5.787e-6;        // |qtop| [m/s]
  const int ncols_src = 30;            // columns 10..39 inside polygontop
  const double area_src = ncols_src * DX * DX;
  const double expect = rate * area_src * 86400.0;
  const double added = store1 - store0;
  std::fprintf(stderr, "  b3: storage added=%.6g m^3, expected recharge=%.6g m^3 (%.2f%%)\n", added,
               expect, 100.0 * added / expect);
  REQUIRE(added > 0.95 * expect);
  REQUIRE(added < 1.05 * expect);

  // (3) A saturated plume formed under the source (head >= 0 somewhere below the recharge band).
  int nsat = 0;
  for (int k = 0; k < NZ; ++k)
    for (int i = 10; i < 40; ++i)
      if (hxz[static_cast<size_t>(k * NX + i)] >= 0.0) ++nsat;
  std::fprintf(stderr, "  b3: saturated cells under source = %d\n", nsat);
  REQUIRE(nsat > 20);

  // (4) Contour match vs the digitized Kirkland reference. Review-tier band: the PCA scheme
  // reproduces the h=0 and h=-400 front positions to within ~0.5 m RMS in a 3 m domain (observed
  // ~0.25-0.29 m). This is NOT a bit-exact legacy/Picard parity claim.
  int n0 = 0, n4 = 0;
  const double rms0 = contourRms(hxz, kRef0, 0.0, n0);
  const double rms4 = contourRms(hxz, kRef400, -400.0, n4);
  std::fprintf(stderr, "  b3: h=0 contour RMS=%.3f m (%d pts), h=-400 RMS=%.3f m (%d pts)\n", rms0,
               n0, rms4, n4);
  REQUIRE(n0 >= 10);
  REQUIRE(n4 >= 8);
  REQUIRE(rms0 < 0.5);
  REQUIRE(rms4 < 0.5);
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
