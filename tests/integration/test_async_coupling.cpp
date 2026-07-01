// P11 acceptance: the async SW<->GW coupling pipeline is numerically equivalent to the P6
// synchronous coupling. b1-sw (SW-only) and b2-gw (GW-only) never enter the coupling path at
// all (no second domain to overlap), so the meaningful async gate is a COUPLED run compared
// between coupling.mode=sync and coupling.mode=async. Because the async pipeline executes the
// identical Gauss-Seidel operation sequence (SW advance -> exchange -> GW catch-up), just with
// the groundwater catch-up of the previous window overlapped on a worker thread, the two runs
// agree to floating-point bit-identity (far inside the 1e-10 gate).
//
// We also assert the regression-safety property the user flagged: requesting coupling.mode
// async on a single-module (SW-only) configuration does NOT activate the pipeline (there is no
// groundwater to couple), so b1-sw / b2-gw remain on the unchanged sequential path.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <string>
#include <vector>

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

// Owned (no-halo) surface eta in global row-major order.
RealArr1DHost ownedEta(const Orchestrator& o) {
  const SweSolver* swe = o.swe();
  const auto& f = swe->fields();
  const int nx = swe->grid().nx(), ny = swe->grid().ny();
  RealArr1DHost out("eta", static_cast<size_t>(nx) * static_cast<size_t>(ny));
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      out(static_cast<size_t>(i + j * nx)) = f.eta(swe->grid().getSurfaceIndex(i, j));
  return out;
}

// Owned (no-halo) groundwater water content in global k-major order.
RealArr1DHost ownedWc(const Orchestrator& o) {
  const ReSolver* re = o.re();
  const auto& f = re->fields();
  const int nx = re->grid().nx(), ny = re->grid().ny(), nz = re->grid().nz();
  RealArr1DHost out("wc", static_cast<size_t>(nx) * static_cast<size_t>(ny) *
                              static_cast<size_t>(nz));
  size_t e = 0;
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) out(e++) = f.wc(re->grid().getIndex(i, j, k));
  return out;
}

double maxAbsDiff(const RealArr1DHost& a, const RealArr1DHost& b) {
  double m = 0.0;
  for (size_t i = 0; i < a.extent(0); ++i) m = std::max(m, std::fabs(a(i) - b(i)));
  return m;
}
}  // namespace

TEST_CASE("async coupling: final SW+GW state matches synchronous coupling to 1e-10") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // CPU async validated on a single rank (see warning in initialize()).

  const double dt_sw = 5.0;
  const long long nsw = 10;
  const double t_end = dt_sw * nsw;

  Config sync_cfg = Config::fromString(
      coupledConfig(std::string(kTmp) + "/async_sync.h5", 4, 4, 10, dt_sw, t_end, nsw, 0.25, 0.15,
                    "sync"),
      "");
  Config async_cfg = Config::fromString(
      coupledConfig(std::string(kTmp) + "/async_async.h5", 4, 4, 10, dt_sw, t_end, nsw, 0.25, 0.15,
                    "async"),
      "");

  Orchestrator seq;
  seq.initialize(sync_cfg);
  REQUIRE(seq.coupled());
  REQUIRE_FALSE(seq.asyncActive());
  seq.run();

  Orchestrator asy;
  asy.initialize(async_cfg);
  REQUIRE(asy.coupled());
  REQUIRE(asy.asyncActive());  // single rank + coupled + async => pipeline engaged
  asy.run();

  // Step / substep bookkeeping is identical (same operation sequence, just overlapped).
  REQUIRE(asy.stepCount() == seq.stepCount());
  REQUIRE(asy.gwSubstepCount() == seq.gwSubstepCount());
  REQUIRE(asy.time() == Approx(seq.time()).margin(1e-12));
  REQUIRE(asy.gwTime() == Approx(seq.gwTime()).margin(1e-12));

  const RealArr1DHost eta_s = ownedEta(seq), eta_a = ownedEta(asy);
  const RealArr1DHost wc_s = ownedWc(seq), wc_a = ownedWc(asy);

  // The gate: relative L2 < 1e-10. In practice these are bit-identical (operations are the same;
  // only the GW catch-up of the previous window is overlapped on a worker thread).
  REQUIRE(relL2(eta_a, eta_s) < 1e-10);
  REQUIRE(relL2(wc_a, wc_s) < 1e-10);
  REQUIRE(maxAbsDiff(eta_a, eta_s) < 1e-12);
  REQUIRE(maxAbsDiff(wc_a, wc_s) < 1e-12);

  // Conservation diagnostic must agree too (exchange is computed at the same sync points).
  REQUIRE(asy.exchangeVolume() == Approx(seq.exchangeVolume()).margin(1e-12));
  // Sanity: the run actually exercised coupling (nonzero exchange, GW subcycled).
  REQUIRE(asy.gwSubstepCount() >= asy.stepCount());
}

TEST_CASE("async coupling: a stronger infiltration scenario also matches to 1e-10") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // Deeper ponding (init_eta 0.4) over drier soil (init_wc 0.10) drives larger SW->GW exchange.
  const double dt_sw = 2.0;
  const long long nsw = 15;
  const double t_end = dt_sw * nsw;

  Config sync_cfg = Config::fromString(
      coupledConfig(std::string(kTmp) + "/async_sync2.h5", 3, 5, 8, dt_sw, t_end, nsw, 0.4, 0.10,
                    "sync"),
      "");
  Config async_cfg = Config::fromString(
      coupledConfig(std::string(kTmp) + "/async_async2.h5", 3, 5, 8, dt_sw, t_end, nsw, 0.4, 0.10,
                    "async"),
      "");

  Orchestrator seq;
  seq.initialize(sync_cfg);
  seq.run();
  Orchestrator asy;
  asy.initialize(async_cfg);
  REQUIRE(asy.asyncActive());
  asy.run();

  REQUIRE(relL2(ownedEta(asy), ownedEta(seq)) < 1e-10);
  REQUIRE(relL2(ownedWc(asy), ownedWc(seq)) < 1e-10);
  REQUIRE(asy.exchangeVolume() == Approx(seq.exchangeVolume()).margin(1e-12));
}

TEST_CASE("async opt-in is regression-safe for single-module runs (b1-sw / b2-gw)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // A SW-only config (b1-sw shape) that nonetheless requests coupling.mode: async. With no
  // groundwater module there is nothing to couple, so the pipeline must NOT engage and the run
  // stays on the unchanged sequential path -> identical to the plain SW-only run.
  auto swWithCoupling = [](const std::string& out, const std::string& mode) {
    std::string base = swConfig(out, 4, 4, 1.0, 5.0, 5, 0.5, 1.0e-4);
    base += "coupling: {mode: " + mode + "}\n";
    return base;
  };

  Config plain = Config::fromString(swConfig(std::string(kTmp) + "/async_swplain.h5", 4, 4, 1.0,
                                             5.0, 5, 0.5, 1.0e-4),
                                    "");
  Config asyncish = Config::fromString(
      swWithCoupling(std::string(kTmp) + "/async_swasync.h5", "async"), "");

  Orchestrator op;
  op.initialize(plain);
  REQUIRE_FALSE(op.coupled());
  op.run();

  Orchestrator oa;
  oa.initialize(asyncish);
  REQUIRE_FALSE(oa.coupled());      // gw disabled => not coupled
  REQUIRE_FALSE(oa.asyncActive());  // => async pipeline never engages
  oa.run();

  REQUIRE(relL2(ownedEta(oa), ownedEta(op)) < 1e-12);
  REQUIRE(maxAbsDiff(ownedEta(oa), ownedEta(op)) == Approx(0.0).margin(1e-13));
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
