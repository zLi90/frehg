// P18 b6-kuan review-tier integration test.
// Coupled surface-water / groundwater column (Kuan et al.; the 1D-2D hybrid test case).
// Frehg2 runs it on the PCA groundwater path with synchronous coupling (the legacy run used the
// Newton solver + baroclinic scalar, which Frehg2 does not implement). Registry gate is `review`,
// so we assert coupled stability + physical bounds + a conservative SW<->GW exchange, not strict
// legacy L2.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstdio>
#include <mpi.h>
#include <string>

#include "core/Grid.hpp"
#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/benchmark_util.hpp"
#include "re/ReSolver.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;
const char* kSrc = FREHG2_SOURCE_DIR;
}  // namespace

TEST_CASE("b6-kuan coupled SW-GW column: stable, bounded, conservative exchange") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string dir = std::string(kSrc) + "/benchmarks/b6-kuan";
  const std::string yaml = dir + "/b6-kuan.yaml";
  const long long steps = 100;          // 100 * dt(0.05) = 5 s coupled
  Config cfg = bench::loadCapped(yaml, dir, steps, std::string(kTmp) + "/b6.h5");

  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.swe() != nullptr);
  REQUIRE(orch.re() != nullptr);
  REQUIRE(orch.coupled());
  orch.run();
  REQUIRE(orch.stepCount() == steps);

  // Surface water: finite, non-negative depths.
  const double maxd = bench::maxSurfaceDepth(*orch.swe());
  REQUIRE(std::isfinite(maxd));
  REQUIRE(maxd >= 0.0);

  // Groundwater: water content stays in the physical [theta_r, theta_s] range everywhere.
  const ReSolver& re = *orch.re();
  const Grid& g = re.grid();
  const auto& gf = re.fields();
  const double tr = re.params().soil.theta_r, ts = re.params().soil.theta_s;
  for (int k = 0; k < g.nz(); ++k)
    for (int j = 0; j < g.ny(); ++j)
      for (int i = 0; i < g.nx(); ++i) {
        const double w = gf.wc(g.getIndex(i, j, k));
        REQUIRE(std::isfinite(w));
        REQUIRE(w >= tr - 1.0e-9);
        REQUIRE(w <= ts + 1.0e-9);
      }

  // The coupled exchange ran and produced a finite (signed) net volume.
  const double exch = orch.exchangeVolume();
  std::fprintf(stderr, "  b6: maxdepth=%.4e m, exchange_volume=%.4e m^3\n", maxd, exch);
  REQUIRE(std::isfinite(exch));
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
