#define FREHG2_TEST_IMPL
#include "frehg2_test.hpp"

#include "io/Config.hpp"
#include "monitoring/MonitorConfig.hpp"

using namespace frehg2;

TEST_CASE("parse monitors.probes and legacy monitoring.points") {
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: true, groundwater: false, solute: false}\n"
      "time: {dt: 1.0, t_end: 10.0}\n"
      "domain: {nx: 5, ny: 5, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1, botz: 0.0}\n"
      "monitors:\n"
      "  probes:\n"
      "    - {name: a, xyz: [2.5, 3.5, 0.0], fields: [eta, u]}\n"
      "  lines:\n"
      "    - {name: gate, p0: [0.0, 2.5, 0.0], p1: [5.0, 2.5, 0.0], field: u}\n"
      "monitoring:\n"
      "  points:\n"
      "    - {i: 1, j: 2, fields: [depth]}\n";
  Config cfg = Config::fromString(yaml, ".");
  const MonitorBundle b = parseMonitors(cfg, true, false);
  REQUIRE(b.probes.size() == 2);
  REQUIRE(b.lines.size() == 1);
  REQUIRE(b.probes[0].name == "a");
  REQUIRE(b.probes[0].fields.size() == 2);
  REQUIRE(b.probes[1].name == "monitor_0");
  REQUIRE(b.lines[0].field == "u");
}

TEST_CASE("parse GW probe with k index") {
  const std::string yaml =
      "schema_version: '2.0'\n"
      "modules: {surface_water: false, groundwater: true, solute: false}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "domain: {nx: 1, ny: 1, nz: 5, dx: 1.0, dy: 1.0, dz: 0.01, botz: -0.05}\n"
      "monitoring:\n"
      "  points:\n"
      "    - {i: 0, j: 0, k: 2, fields: [head]}\n";
  Config cfg = Config::fromString(yaml, ".");
  const MonitorBundle b = parseMonitors(cfg, false, true);
  REQUIRE(b.probes.size() == 1);
  REQUIRE(b.probes[0].subsurface == true);
  REQUIRE(b.probes[0].gk == 2);
}
