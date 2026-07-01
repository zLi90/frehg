// P18 b4-govindaraju review-tier integration test.
// The legacy SERGHEI case is a rainfall-driven overland-flow plane with Chezy friction; Frehg2
// runs it with Manning friction + a polygon critical-depth outlet. The registry gates b4
// `review` (tolerance: null), so we assert review-tier physics: the run is stable, depths stay
// finite and non-negative, rainfall wets the plane, and no water is spuriously created (storage
// + outlet outflow <= rainfall delivered). This is NOT a strict legacy-L2 reproduction.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cstdio>
#include <mpi.h>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/benchmark_util.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;
const char* kSrc = FREHG2_SOURCE_DIR;
}  // namespace

TEST_CASE("b4-govindaraju overland flow: stable, mass-conserving, wets the plane") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // review-tier benchmark runs serially

  const std::string dir = std::string(kSrc) + "/benchmarks/b4-govindaraju";
  const std::string yaml = dir + "/b4-govindaraju.yaml";
  const long long steps = 600;          // 600 * dt(0.05) = 30 s of rainfall
  Config cfg = bench::loadCapped(yaml, dir, steps, std::string(kTmp) + "/b4.h5");

  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.swe() != nullptr);
  orch.run();
  REQUIRE(orch.stepCount() == steps);

  const double vol = bench::surfaceWaterVolume(*orch.swe());
  const double maxd = bench::maxSurfaceDepth(*orch.swe());
  const double outflow = orch.polygonOutflowVolume();

  // Finite, non-negative; rainfall entered the system (became storage and/or left via the
  // outlet — on this steep, finely-celled plane most of the early water drains downslope).
  REQUIRE(std::isfinite(vol));
  REQUIRE(std::isfinite(maxd));
  REQUIRE(maxd >= 0.0);
  REQUIRE(vol >= 0.0);
  REQUIRE(outflow >= 0.0);
  REQUIRE(vol + outflow > 0.0);

  // No spurious creation: storage + outlet loss <= rainfall delivered over the full plane.
  const double dx = 0.109725, dy = 0.109725;
  const double area = 200.0 * 10.0 * dx * dy;
  const double max_rate = 2.8222194e-05;       // peak hyetograph value [m/s]
  const double rain_upper = max_rate * (steps * 0.05) * area;
  std::fprintf(stderr,
               "  b4: storage=%.4e m^3, outflow=%.4e m^3, rain_upper=%.4e m^3, maxdepth=%.3e m\n",
               vol, outflow, rain_upper, maxd);
  REQUIRE(vol + outflow <= rain_upper * 1.001);
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
