// P17 acceptance: v1/experimental -> frozen v2 YAML migration.
// Loads the archived legacy/benchmarks/*/input.v1.yaml fixtures, migrates them, and asserts the
// frozen v2 keys are produced. Also checks idempotency on a real v2 benchmark fixture and the
// boundary/source `type` rename. (No Orchestrator here: this links Frehg2::io only.)
#define FREHG2_TEST_IMPL
#include <string>

#include "frehg2_test.hpp"
#include "io/Config.hpp"
#include "io/YamlMigration.hpp"

using namespace frehg2;

namespace {
const char* kSrcDir = FREHG2_SOURCE_DIR;
const char* kBenchDir = FREHG2_BENCH_DIR;
std::string srcPath(const std::string& rel) { return std::string(kSrcDir) + "/" + rel; }
std::string benchPath(const std::string& rel) { return std::string(kBenchDir) + "/" + rel; }

// Migrate a v1 fixture file and return it parsed as a (validated) Config.
Config migrateFixture(const std::string& v1_rel) {
  YAML::Node v1 = YAML::LoadFile(srcPath(v1_rel));
  YAML::Node v2 = migrateV1ToV2(v1);
  std::stringstream ss;
  ss << v2;
  return Config::fromString(ss.str());
}
}  // namespace

TEST_CASE("migrate b1-sw v1 -> v2: grid/time/io renames and required fields") {
  Config c = migrateFixture("legacy/benchmarks/b1-sw/input.v1.yaml");

  REQUIRE(c.get<std::string>("schema_version") == std::string("2.0"));
  REQUIRE(c.has("simulation"));
  REQUIRE(c.has("modules"));
  REQUIRE(c.has("output"));

  // grid.* -> domain.*  (+ bot_z -> botz)
  REQUIRE(c.get<int>("domain.nx") == 1);
  REQUIRE(c.get<int>("domain.ny") == 10);
  REQUIRE(c.get<int>("domain.nz") == 30);
  REQUIRE(c.get<double>("domain.botz") == Approx(-3.0).margin(1e-12));
  REQUIRE(c.get<double>("domain.dx") == Approx(80.0).margin(1e-12));
  REQUIRE_FALSE(c.has("grid"));
  REQUIRE_FALSE(c.has("domain.bot_z"));

  // time.{Tend,max_step,dt_out} -> frozen names
  REQUIRE(c.get<double>("time.t_end") == Approx(18000.0).margin(1e-6));
  REQUIRE(c.get<int>("time.max_steps") == 100);
  REQUIRE(c.get<double>("time.output_interval") == Approx(1800.0).margin(1e-9));
  REQUIRE_FALSE(c.has("time.Tend"));
  REQUIRE_FALSE(c.has("time.dt_out"));

  // io.dir -> output.filename
  REQUIRE(c.get<std::string>("output.format") == std::string("hdf5"));
  REQUIRE(c.get<std::string>("output.filename") == std::string("out/b1-sw.h5"));
  REQUIRE_FALSE(c.has("io"));

  REQUIRE(c.get<bool>("modules.surface_water") == true);
}

TEST_CASE("migrate b2-gw v1 -> v2: soil.uniform -> frozen soil.types[]") {
  Config c = migrateFixture("legacy/benchmarks/b2-gw/input.v1.yaml");

  REQUIRE(c.get<std::string>("schema_version") == std::string("2.0"));
  REQUIRE(c.get<int>("domain.nz") == 100);
  REQUIRE(c.get<double>("domain.botz") == Approx(-1.0).margin(1e-12));
  REQUIRE(c.get<double>("time.t_end") == Approx(46800.0).margin(1e-6));
  REQUIRE(c.get<double>("time.output_interval") == Approx(11700.0).margin(1e-6));

  // soil.uniform (+ flat params) -> a single explicit class. Sequence entries are read from
  // the raw node (the dotted Config path does not index sequences).
  REQUIRE_FALSE(c.has("soil.uniform"));
  REQUIRE(c.indexedCount("soil.types") == 1u);
  const YAML::Node cls = c.root()["soil"]["types"][0];
  REQUIRE(cls["theta_s"].as<double>() == Approx(0.33).margin(1e-12));
  REQUIRE(cls["theta_r"].as<double>() == Approx(0.0).margin(1e-12));
  REQUIRE(cls["vg"]["alpha"].as<double>() == Approx(1.43).margin(1e-12));
  REQUIRE(cls["vg"]["n"].as<double>() == Approx(1.56).margin(1e-12));
  REQUIRE(cls["k_sat"]["z"].as<double>() == Approx(2.89e-6).epsilon(1e-9));
  REQUIRE(c.get<bool>("modules.groundwater") == true);
}

TEST_CASE("migration is idempotent on a production v2 fixture") {
  YAML::Node v2 = YAML::LoadFile(benchPath("b2-gw/b2-gw.yaml"));
  YAML::Node again = migrateV1ToV2(v2);
  std::stringstream ss;
  ss << again;
  Config c = Config::fromString(ss.str());

  // Frozen keys untouched; no v1 artefacts introduced.
  REQUIRE(c.get<std::string>("schema_version") == std::string("2.0"));
  REQUIRE(c.get<int>("domain.nz") == 100);
  REQUIRE(c.get<double>("time.t_end") == Approx(46800.0).margin(1e-6));
  REQUIRE(c.indexedCount("soil.types") == v2["soil"]["types"].size());
  REQUIRE_FALSE(c.has("grid"));
  REQUIRE_FALSE(c.has("io"));
  REQUIRE_FALSE(c.has("time.Tend"));
}

TEST_CASE("boundary/source `type` kinds are renamed to frozen vocabulary") {
  const std::string v1 =
      "simulation: {id: bc_test, mode: surface_water}\n"
      "domain: {nx: 1, ny: 1, nz: 1}\n"
      "time: {dt: 1.0, Tend: 1.0}\n"
      "modules: {surface_water: true}\n"
      "boundaries:\n"
      "  - {name: outlet, type: discharge, rate: 1.0, vertices: [[0,0],[1,0],[1,1]]}\n"
      "  - {name: weir, type: critical, vertices: [[0,0],[1,0],[1,1]]}\n"
      "sources:\n"
      "  - {name: src, type: inflow, rate: 2.0, vertices: [[0,0],[1,0],[1,1]]}\n"
      "output: {format: hdf5, directory: results}\n";
  Config c = Config::fromString(migrateV1ToV2String(v1));

  REQUIRE(c.root()["boundaries"][0]["type"].as<std::string>() == std::string("bc_discharge"));
  REQUIRE(c.root()["boundaries"][1]["type"].as<std::string>() == std::string("bc_critical"));
  REQUIRE(c.root()["sources"][0]["type"].as<std::string>() == std::string("inflow_rate"));
  // output.directory -> output.filename (joined with sim id).
  REQUIRE(c.get<std::string>("output.filename") == std::string("results/bc_test.h5"));
  REQUIRE(c.get<double>("time.t_end") == Approx(1.0).margin(1e-12));
}
