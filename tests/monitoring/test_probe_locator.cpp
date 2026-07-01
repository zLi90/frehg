#define FREHG2_TEST_IMPL
#include "frehg2_test.hpp"

#include "core/Grid.hpp"
#include "monitoring/ProbeLocator.hpp"

using namespace frehg2;

TEST_CASE("nearest-cell lookup from xyz on structured grid") {
  MonitorBundle bundle;
  ProbeSpec p;
  p.name = "center";
  p.coords_from_xyz = true;
  p.x = 45.0;
  p.y = 365.0;
  p.z = 0.0;
  bundle.probes.push_back(p);

  ProbeSpec legacy;
  legacy.name = "legacy";
  legacy.indices_from_grid = true;
  legacy.gi = 0;
  legacy.gj = 4;
  bundle.probes.push_back(legacy);

  buildProbeLocations(bundle, 1, 10, 1, 80.0, 80.0, 0.1, 0.0, 0.0, -3.0, nullptr);

  REQUIRE(bundle.probes[0].gi == 0);
  REQUIRE(bundle.probes[0].gj == 4);
  REQUIRE(bundle.probes[1].gi == 0);
  REQUIRE(bundle.probes[1].gj == 4);
}

TEST_CASE("xyz maps to expected cell on uniform 10x10 grid") {
  MonitorBundle bundle;
  ProbeSpec p;
  p.name = "p";
  p.coords_from_xyz = true;
  p.x = 45.0;
  p.y = 25.0;
  bundle.probes.push_back(p);
  buildProbeLocations(bundle, 10, 10, 1, 10.0, 10.0, 1.0, 0.0, 0.0, 0.0, nullptr);
  REQUIRE(bundle.probes[0].gi == 4);
  REQUIRE(bundle.probes[0].gj == 2);
}
