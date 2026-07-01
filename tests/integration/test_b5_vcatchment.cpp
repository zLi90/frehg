// P18 b5-vcatchment review-tier integration test.
// V-catchment overland flow on the real SERGHEI DEM (101 x 55), Manning friction, polygon outlet.
// Registry gate is `review` (the legacy refs are an inter-model CATHY/HGS/ParFlow comparison), so
// we assert review-tier physics: stable run, finite non-negative depths, rainfall wets the
// catchment, and no spurious water creation. SW-only (subsurface leg deferred; see the doc).
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

TEST_CASE("b5-vcatchment overland flow: stable, mass-conserving on the real DEM") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string dir = std::string(kSrc) + "/benchmarks/b5-vcatchment";
  const std::string yaml = dir + "/b5-vcatchment.yaml";
  const long long steps = 200;          // 200 * dt(1.0) = 200 s of rainfall
  Config cfg = bench::loadCapped(yaml, dir, steps, std::string(kTmp) + "/b5.h5");

  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.swe() != nullptr);
  REQUIRE(orch.swe()->grid().nx() == 101);
  REQUIRE(orch.swe()->grid().ny() == 55);
  orch.run();
  REQUIRE(orch.stepCount() == steps);

  const double vol = bench::surfaceWaterVolume(*orch.swe());
  const double maxd = bench::maxSurfaceDepth(*orch.swe());
  const double outflow = orch.polygonOutflowVolume();

  REQUIRE(std::isfinite(vol));
  REQUIRE(std::isfinite(maxd));
  REQUIRE(maxd >= 0.0);
  REQUIRE(vol >= 0.0);
  REQUIRE(outflow >= 0.0);
  REQUIRE(vol + outflow > 0.0);

  const double area = 101.0 * 55.0 * 1.0 * 1.0;
  const double rate = 2.7777778e-05;            // 100 mm/h [m/s]
  const double rain_upper = rate * (steps * 1.0) * area;
  std::fprintf(stderr, "  b5: storage=%.4e m^3, outflow=%.4e m^3, rain_upper=%.4e m^3\n", vol,
               outflow, rain_upper);
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
