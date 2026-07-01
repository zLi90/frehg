// P17 acceptance: the production driver enforces the frozen v2 schema.
//   * Orchestrator::initialize() refuses a missing/incorrect schema_version.
//   * Orchestrator::initialize() refuses a missing required top-level section.
//   * A v1 fixture migrated by migrateV1ToV2 runs end-to-end on the production driver
//     (proving the v1 -> v2 round-trip is actually runnable, not just parseable).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <sstream>
#include <string>

#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "io/YamlMigration.hpp"
#include "re/ReSolver.hpp"

using namespace frehg2;
using namespace frehg2::orch_test;

namespace {
const char* kTmp = FREHG2_IO_TMP;
const char* kSrcDir = FREHG2_SOURCE_DIR;

// Remove the first line that contains `token` from a YAML string.
std::string dropLine(const std::string& yaml, const std::string& token) {
  std::istringstream in(yaml);
  std::ostringstream out;
  std::string line;
  bool dropped = false;
  while (std::getline(in, line)) {
    if (!dropped && line.find(token) != std::string::npos) {
      dropped = true;
      continue;
    }
    out << line << "\n";
  }
  return out.str();
}

bool initThrows(const std::string& yaml, const std::string& expect_substr) {
  bool threw = false;
  try {
    Config cfg = Config::fromString(yaml, "");
    Orchestrator orch;
    orch.initialize(cfg);
  } catch (const std::runtime_error& e) {
    threw = true;
    REQUIRE(std::string(e.what()).find(expect_substr) != std::string::npos);
  }
  return threw;
}
}  // namespace

TEST_CASE("Orchestrator refuses missing/incorrect schema_version") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string good = swConfig(std::string(kTmp) + "/schema_ok.h5", 3, 3, 1.0, 1.0, 1, 0.5);

  // Good config initializes.
  {
    Config cfg = Config::fromString(good, "");
    Orchestrator orch;
    orch.initialize(cfg);
    REQUIRE(orch.swe() != nullptr);
  }

  // Missing schema_version.
  REQUIRE(initThrows(dropLine(good, "schema_version"), "schema_version"));

  // Wrong schema_version.
  std::string wrong = good;
  const auto pos = wrong.find("'2.0'");
  REQUIRE(pos != std::string::npos);
  wrong.replace(pos, 5, "'1.0'");
  REQUIRE(initThrows(wrong, "schema_version"));
}

TEST_CASE("Orchestrator refuses a missing required top-level section") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string good = swConfig(std::string(kTmp) + "/schema_sect.h5", 3, 3, 1.0, 1.0, 1, 0.5);
  // 'output:' is required by the production gate.
  REQUIRE(initThrows(dropLine(good, "output:"), "output"));
}

TEST_CASE("migrated v1 fixture runs on the production driver (b2-gw)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const std::string v1_path =
      std::string(kSrcDir) + "/legacy/benchmarks/b2-gw/input.v1.yaml";
  YAML::Node v1 = YAML::LoadFile(v1_path);
  YAML::Node v2 = migrateV1ToV2(v1);
  // Cap the run so the test is fast; the migrated config is otherwise complete and runnable.
  v2["time"]["max_steps"] = 20;
  std::stringstream ss;
  ss << v2;

  Config cfg = Config::fromString(ss.str(), kTmp);
  Orchestrator orch;
  orch.initialize(cfg);
  REQUIRE(orch.re() != nullptr);
  REQUIRE(orch.re()->grid().nz() == 100);
  orch.run();
  REQUIRE(orch.stepCount() == 20);
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
