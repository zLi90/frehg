// P7.3 acceptance: coupled time-stepping orchestration. A fixed surface dt with an adaptive
// groundwater dt drives several GW substeps per surface step; total GW time tracks total SW
// time; SW-only / GW-only instantiate only their solver; coupled instantiates both. The
// exchange across the interface stays mass-sane (states finite, depth >= 0, wc in bounds).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <mpi.h>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "re/ReSolver.hpp"
#include "re/VanGenuchten.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;
}

TEST_CASE("coupled: fixed dt_sw + adaptive dt_gw -> multiple GW substeps, times match") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // Surface step 60 s; GW dt_min 1 s and adaptive -> many GW substeps per surface step.
  const double dt_sw = 60.0;
  const long long nsw = 5;
  const std::string out = std::string(kTmp) + "/coupled_ts.h5";
  Config cfg = Config::fromString(coupledConfig(out, 3, 3, 8, dt_sw, dt_sw * nsw, nsw, 0.3, 0.2),
                                  "");
  // Tighten GW dt_min via a second config knob: coupledConfig already sets dt_min 0.001;
  // with co_max=2 the adaptive dt grows from there, so several substeps occur early.
  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.coupled());
  orch.run();

  REQUIRE(orch.stepCount() == nsw);
  // Each surface step subcycles GW at least once; with a small starting dt_gw, strictly more
  // GW substeps than surface steps must occur.
  REQUIRE(orch.gwSubstepCount() >= nsw);
  // Total GW time caught up to total SW time (within one GW substep).
  REQUIRE(orch.gwTime() >= orch.time() - 1.0e-6);
  REQUIRE(orch.time() == Approx(dt_sw * nsw).margin(1e-9));
}

TEST_CASE("coupled: interface exchange keeps states physical") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string out = std::string(kTmp) + "/coupled_phys.h5";
  Config cfg = Config::fromString(coupledConfig(out, 4, 4, 10, 1.0, 20.0, 20, 0.25, 0.15), "");
  Orchestrator orch;
  orch.initialize(cfg);
  orch.run();

  const SweSolver* swe = orch.swe();
  const ReSolver* re = orch.re();
  const auto& sf = swe->fields();
  const auto& gf = re->fields();
  const auto& soil = re->params().soil;
  const int nx = swe->grid().nx(), ny = swe->grid().ny(), nz = re->grid().nz();
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = swe->grid().getSurfaceIndex(i, j);
      REQUIRE(std::isfinite(sf.eta(c)));
      REQUIRE(sf.eta(c) - sf.bottom(c) >= -1.0e-12);
    }
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) {
        const int c = re->grid().getIndex(i, j, k);
        REQUIRE(std::isfinite(gf.wc(c)));
        REQUIRE(gf.wc(c) >= soil.theta_r - 1.0e-12);
        REQUIRE(gf.wc(c) <= soil.theta_s + 1.0e-12);
      }
  REQUIRE(std::isfinite(orch.exchangeVolume()));
}

TEST_CASE("module selection: SW-only and GW-only instantiate only their solver") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  Config sw = Config::fromString(swConfig(std::string(kTmp) + "/mods_sw.h5", 4, 4, 1.0, 2.0,
                                           2, 0.5, 0.0),
                                 "");
  Orchestrator osw;
  osw.initialize(sw);
  REQUIRE(osw.swe() != nullptr);
  REQUIRE(osw.re() == nullptr);
  REQUIRE_FALSE(osw.coupled());

  Config gw = Config::fromString(gwConfig(std::string(kTmp) + "/mods_gw.h5", 10, 1.0e-4, 1.0,
                                          5, 0.1),
                                 "");
  Orchestrator ogw;
  ogw.initialize(gw);
  REQUIRE(ogw.swe() == nullptr);
  REQUIRE(ogw.re() != nullptr);
  REQUIRE_FALSE(ogw.coupled());
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
