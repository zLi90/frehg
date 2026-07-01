// P7.4 acceptance: checkpoint + restart. Running to completion in one shot vs. running to a
// mid-point, checkpointing, then restarting and continuing must agree to L2 < 1e-10 (same
// executable/backend). Restarting from a non-existent file gives a clear error.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <sstream>
#include <string>

#include "core/Grid.hpp"
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
const std::string kB1Dir = std::string(FREHG2_SOURCE_DIR) + "/benchmarks/b1-sw";

RealArr1DHost ownedDepth(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const int nx = g.nx(), ny = g.ny();
  RealArr1DHost d("d", static_cast<size_t>(nx * ny));
  const auto& f = swe.fields();
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      const int c = g.getSurfaceIndex(i, j);
      d(static_cast<size_t>(i + j * nx)) = std::max(f.eta(c) - f.bottom(c), 0.0);
    }
  return d;
}

// Mirror of Hdf5Writer::timeKey so tests can predict checkpoint file names: whole seconds
// stay plain integers; sub-second times get a lossless trimmed decimal.
std::string ckptKey(double t) {
  const long long r = std::llround(t);
  if (std::fabs(t - static_cast<double>(r)) < 1.0e-9) return std::to_string(r);
  std::ostringstream os;
  os << std::setprecision(12) << t;
  return os.str();
}

RealArr1DHost wcColumn(const ReSolver& re) {
  const Grid& g = re.grid();
  const int nz = g.nz();
  RealArr1DHost w("w", static_cast<size_t>(nz));
  for (int k = 0; k < nz; ++k) w(static_cast<size_t>(k)) = re.fields().wc(g.getIndex(0, 0, k));
  return w;
}
}  // namespace

TEST_CASE("b1-sw restart matches continuous run (L2 < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const double dt = 5.0;
  const double t_mid = 50.0;   // 10 steps
  const double t_end = 100.0;  // 20 steps total
  const double init_eta = 0.7;  // wet basin -> real flow + drag

  // Continuous reference.
  {
    const std::string out = std::string(kTmp) + "/restart_cont.h5";
    Config c = Config::fromString(b1Config(out, init_eta, dt, t_end, 20), kB1Dir);
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
  }
  RealArr1DHost cont;
  {
    const std::string out = std::string(kTmp) + "/restart_cont2.h5";
    Config c = Config::fromString(b1Config(out, init_eta, dt, t_end, 20), kB1Dir);
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    cont = ownedDepth(*orch.swe());
  }

  // Phase A: run to the mid-point and checkpoint there.
  const std::string outA = std::string(kTmp) + "/restart_phaseA.h5";
  {
    Config c = Config::fromString(b1Config(outA, init_eta, dt, t_mid, 10, t_mid), kB1Dir);
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
  }
  const std::string ckpt = std::string(kTmp) + "/restart_phaseA.ckpt.50.h5";

  // Phase B: restart from the checkpoint and continue to the end.
  RealArr1DHost restarted;
  {
    Config c = Config::fromString(b1Config(outA, init_eta, dt, t_end, 20), kB1Dir);
    Orchestrator orch;
    orch.initialize(c);
    orch.restart(ckpt, t_mid);
    restarted = ownedDepth(*orch.swe());
    REQUIRE(orch.time() == Approx(t_end).margin(1e-9));
  }

  const double l2 = relL2(restarted, cont);
  std::fprintf(stderr, "  b1-sw restart-vs-continuous rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-10);
}

// Build a GW config with a dt_checkpoint injected into the time map.
namespace {
std::string gwConfigWithCheckpoint(const std::string& out, int nz, double dt, double t_end,
                                   long long max_steps, double init_wc, double dt_checkpoint) {
  std::string yaml = gwConfig(out, nz, dt, t_end, max_steps, init_wc);
  const std::string needle = "output_interval: 0}";
  const std::string repl = "output_interval: 0, dt_checkpoint: " +
                           std::to_string(dt_checkpoint) + ", max_checkpoints: 4}";
  const auto pos = yaml.find(needle);
  if (pos != std::string::npos) yaml.replace(pos, needle.size(), repl);
  return yaml;
}
}  // namespace

TEST_CASE("b2-gw restart matches continuous run (L2 < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nz = 40;
  const double dt = 1.0e-4;
  const double t_end = 1.0e6;  // governed by max_steps

  // Continuous reference: 60 adaptive steps.
  RealArr1DHost cont;
  {
    Config c = Config::fromString(gwConfig(std::string(kTmp) + "/rgw_cont.h5", nz, dt, t_end,
                                           60, 0.033),
                                  "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    cont = wcColumn(*orch.re());
  }

  // Probe: learn the adaptive simulated time after 30 steps (deterministic).
  double t30 = 0.0;
  {
    Config c = Config::fromString(gwConfig(std::string(kTmp) + "/rgw_probe.h5", nz, dt, t_end,
                                           30, 0.033),
                                  "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
    t30 = orch.time();
  }

  // Phase A: run 30 steps and checkpoint exactly at t30.
  const std::string outA = std::string(kTmp) + "/rgw_phaseA.h5";
  {
    Config c = Config::fromString(gwConfigWithCheckpoint(outA, nz, dt, t_end, 30, 0.033, t30),
                                  "");
    Orchestrator orch;
    orch.initialize(c);
    orch.run();
  }
  const std::string ckpt = std::string(kTmp) + "/rgw_phaseA.ckpt." + ckptKey(t30) + ".h5";

  // Phase B: restart from the checkpoint and run the remaining 30 steps (60 total). The GW
  // run length is governed by max_steps, and phase B's step counter starts at zero, so it is
  // set to the remaining step count.
  RealArr1DHost restarted;
  {
    Config c = Config::fromString(gwConfig(outA, nz, dt, t_end, 30, 0.033), "");
    Orchestrator orch;
    orch.initialize(c);
    orch.restart(ckpt, t30);
    restarted = wcColumn(*orch.re());
  }

  const double l2 = relL2(restarted, cont);
  std::fprintf(stderr, "  b2-gw restart-vs-continuous rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-10);
}

TEST_CASE("restart from a non-existent checkpoint gives a clear error") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  Config c = Config::fromString(swConfig(std::string(kTmp) + "/restart_err.h5", 4, 4, 1.0, 5.0,
                                         5, 0.5, 0.0),
                                "");
  Orchestrator orch;
  orch.initialize(c);
  bool threw = false;
  std::string msg;
  try {
    orch.restart(std::string(kTmp) + "/does_not_exist.ckpt.0.h5", 0.0);
  } catch (const std::runtime_error& e) {
    threw = true;
    msg = e.what();
  }
  REQUIRE(threw);
  REQUIRE(msg.find("does_not_exist") != std::string::npos);
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
