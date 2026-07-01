// P3.1 acceptance: YAML Config parsing, typed access, schema validation, path resolution.
#define FREHG2_TEST_IMPL
#include <string>

#include "frehg2_test.hpp"
#include "io/Config.hpp"

using namespace frehg2;

namespace {
const char* kBenchDir = FREHG2_BENCH_DIR;  // set by CMake to the benchmarks/ directory
std::string benchPath(const std::string& rel) { return std::string(kBenchDir) + "/" + rel; }
}  // namespace

TEST_CASE("parse b1-sw.yaml: required fields and types") {
  Config cfg = Config::fromFile(benchPath("b1-sw/b1-sw.yaml"));
  REQUIRE(cfg.get<std::string>("simulation.id") == std::string("b1-sw"));
  REQUIRE(cfg.get<int>("domain.nx") == 1);
  REQUIRE(cfg.get<int>("domain.ny") == 10);
  REQUIRE(cfg.get<int>("domain.nz") == 30);
  REQUIRE(cfg.get<double>("time.dt") == Approx(5.0).margin(1e-12));
  REQUIRE(cfg.get<double>("time.t_end") == Approx(18000.0).margin(1e-9));
  REQUIRE(cfg.get<bool>("modules.surface_water") == true);
  REQUIRE(cfg.get<bool>("modules.groundwater") == false);
  REQUIRE(cfg.get<double>("surface_water.gravity") == Approx(9.81).margin(1e-12));
}

TEST_CASE("parse b2-gw.yaml: required fields") {
  Config cfg = Config::fromFile(benchPath("b2-gw/b2-gw.yaml"));
  REQUIRE(cfg.get<int>("domain.nx") > 0);
  REQUIRE(cfg.get<int>("domain.nz") > 0);
  REQUIRE(cfg.has("time.t_end"));
  REQUIRE(cfg.has("groundwater.solver"));
}

TEST_CASE("missing required key throws with key name") {
  Config cfg = Config::fromFile(benchPath("b1-sw/b1-sw.yaml"));
  bool threw = false;
  try {
    (void)cfg.get<int>("domain.does_not_exist");
  } catch (const std::runtime_error& e) {
    threw = true;
    const std::string msg = e.what();
    REQUIRE(msg.find("domain.does_not_exist") != std::string::npos);
  }
  REQUIRE(threw);
}

TEST_CASE("getOr returns default for missing, value for present") {
  Config cfg = Config::fromFile(benchPath("b1-sw/b1-sw.yaml"));
  REQUIRE(cfg.getOr<int>("domain.not_here", 42) == 42);
  REQUIRE(cfg.getOr<int>("domain.nx", 999) == 1);
  REQUIRE(cfg.getOr<double>("surface_water.manning", -1.0) == Approx(0.019).margin(1e-9));
}

TEST_CASE("wrong type throws") {
  Config cfg = Config::fromFile(benchPath("b1-sw/b1-sw.yaml"));
  bool threw = false;
  try {
    (void)cfg.get<int>("simulation.id");  // string -> int
  } catch (const std::runtime_error&) {
    threw = true;
  }
  REQUIRE(threw);
}

TEST_CASE("indexedCount on sequences") {
  Config cfg = Config::fromFile(benchPath("b1-sw/b1-sw.yaml"));
  REQUIRE(cfg.indexedCount("groundwater.bc_type_gw") == 6u);
  REQUIRE(cfg.indexedCount("soil.types") == 1u);
  REQUIRE(cfg.indexedCount("does.not.exist") == 0u);
}

TEST_CASE("relative path resolved against config directory") {
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: true}\n"
      "time: {dt: 1.0, t_end: 10.0}\n"
      "domain: {nx: 2, ny: 3, nz: 4}\n"
      "output: {bathymetry: dem/bath.asc}\n";
  Config cfg = Config::fromString(yaml, "/data/run1");
  REQUIRE(cfg.resolvePath(cfg.get<std::string>("output.bathymetry")) ==
          std::string("/data/run1/dem/bath.asc"));
  REQUIRE(cfg.resolvePath("/abs/path.asc") == std::string("/abs/path.asc"));
}

TEST_CASE("schema validation rejects missing/invalid required sections") {
  bool threw = false;
  try {
    Config::fromString("time: {dt: 1.0, t_end: 2.0}\ndomain: {nx: 1, ny: 1, nz: 1}\n");
  } catch (const std::runtime_error&) {
    threw = true;  // 'modules' missing
  }
  REQUIRE(threw);

  threw = false;
  try {
    Config::fromString(
        "modules: {surface_water: true}\ntime: {dt: 1.0, t_end: 2.0}\n"
        "domain: {nx: -1, ny: 1, nz: 1}\n");
  } catch (const std::runtime_error&) {
    threw = true;  // nx not positive
  }
  REQUIRE(threw);
}
