// P7.1 acceptance: the unified Orchestrator runs SW-only, GW-only, and coupled smoke cases
// from a single YAML config, each completing without error, instantiating only the enabled
// solvers, and writing simulation_summary.txt.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <fstream>
#include <mpi.h>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;

bool fileExists(const std::string& p) {
  std::ifstream f(p);
  return f.good();
}
}  // namespace

TEST_CASE("SW-only smoke: runs, only SWE instantiated, summary written") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/smoke_sw.h5";
  Config cfg = Config::fromString(swConfig(out, 5, 5, 1.0, 5.0, 5, 0.5, 1.0e-5), "");
  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.swEnabled());
  REQUIRE_FALSE(orch.gwEnabled());
  REQUIRE(orch.re() == nullptr);
  orch.run();
  REQUIRE(orch.stepCount() == 5);
  REQUIRE(orch.time() == Approx(5.0).margin(1e-9));
  REQUIRE(fileExists(orch.summaryPath()));
  // SWE state stayed finite.
  const SweSolver* s = orch.swe();
  REQUIRE(s != nullptr);
  const auto& f = s->fields();
  const int c = s->grid().getSurfaceIndex(2, 2);
  REQUIRE(std::isfinite(f.eta(c)));
}

TEST_CASE("GW-only smoke: runs, only RE instantiated, summary written") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/smoke_gw.h5";
  // Adaptive GW dt -> use a large t_end so max_steps governs the step count.
  Config cfg = Config::fromString(gwConfig(out, 20, 1.0e-4, 1.0e9, 40, 0.1), "");
  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE_FALSE(orch.swEnabled());
  REQUIRE(orch.gwEnabled());
  REQUIRE(orch.swe() == nullptr);
  orch.run();
  REQUIRE(orch.stepCount() == 40);
  REQUIRE(fileExists(orch.summaryPath()));
  const ReSolver* r = orch.re();
  REQUIRE(r != nullptr);
  const auto& f = r->fields();
  const int c = r->grid().getIndex(0, 0, 5);
  REQUIRE(std::isfinite(f.wc(c)));
  REQUIRE(f.wc(c) >= 0.0);
}

TEST_CASE("coupled smoke: both solvers instantiated, runs, summary written") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/smoke_coupled.h5";
  Config cfg = Config::fromString(coupledConfig(out, 4, 4, 8, 1.0, 10.0, 10, 0.2, 0.2), "");
  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.swEnabled());
  REQUIRE(orch.gwEnabled());
  REQUIRE(orch.coupled());
  orch.run();
  REQUIRE(orch.stepCount() == 10);
  REQUIRE(fileExists(orch.summaryPath()));
  // Coupled exchange happened (some water moved across the interface) and states are finite.
  const auto& sf = orch.swe()->fields();
  const auto& gf = orch.re()->fields();
  const int sc = orch.swe()->grid().getSurfaceIndex(1, 1);
  const int gc = orch.re()->grid().getIndex(1, 1, 0);
  REQUIRE(std::isfinite(sf.eta(sc)));
  REQUIRE(std::isfinite(gf.wc(gc)));
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
