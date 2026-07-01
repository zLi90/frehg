// P12.3.5 acceptance: the YAML polygon schema parses 5+ polygon configurations and dispatches to
// the correct BC / source kinds, vertices, and rates. Also asserts fail-loud behavior on a bad
// type and too-few vertices, and that the legacy scalar `sources:` map is left untouched.
#define FREHG2_TEST_IMPL
#include <string>

#include "bc/PolygonConfig.hpp"
#include "frehg2_test.hpp"
#include "io/Config.hpp"

using namespace frehg2;

namespace {
const char* kCfg = R"(schema_version: '2.0'
simulation: {id: poly, mode: surface_water}
domain: {nx: 10, ny: 10, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: true, groundwater: true, solute: false}
boundaries:
  - {name: outlet, type: bc_discharge, vertices: [[9,0],[10,0],[10,10],[9,10]], rate: 2.5}
  - {name: wall,   type: bc_depth,     vertices: [[0,0],[1,0],[1,10],[0,10]], rate: 0.3}
  - {name: weir,   type: bc_critical,  vertices: [[4,9],[6,9],[6,10],[4,10]]}
sources:
  - {name: in,   type: inflow_rate,     vertices: [[0,0],[2,0],[2,2],[0,2]], rate: 1.0}
  - {name: rainpatch, type: rainfall_rate, vertices: [[5,5],[7,5],[7,7],[5,7]], rate: 1.0e-5}
  - {name: well, type: extraction_well, vertices: [[8,8],[9,8],[9,9],[8,9]], rate: 0.05}
)";
}  // namespace

TEST_CASE("yaml polygon: parses boundaries and sources with correct dispatch") {
  Config cfg = Config::fromString(kCfg, "");
  PolygonRegions pr = parsePolygonRegions(cfg);

  REQUIRE(pr.boundaries.size() == 3u);
  REQUIRE(pr.sources.size() == 3u);  // total 6 polygons (>= 5 required)

  REQUIRE(pr.boundaries[0].kind == BcKind::Discharge);
  REQUIRE(pr.boundaries[0].polygon.name == "outlet");
  REQUIRE(pr.boundaries[0].rate == Approx(2.5).margin(1e-12));
  REQUIRE(pr.boundaries[0].polygon.vertices.size() == 4u);
  REQUIRE(pr.boundaries[1].kind == BcKind::Depth);
  REQUIRE(pr.boundaries[1].rate == Approx(0.3).margin(1e-12));
  REQUIRE(pr.boundaries[2].kind == BcKind::Critical);
  REQUIRE(pr.boundaries[2].rate == Approx(0.0).margin(1e-12));  // default when absent

  REQUIRE(pr.sources[0].kind == SourceKind::InflowRate);
  REQUIRE(pr.sources[0].rate == Approx(1.0).margin(1e-12));
  REQUIRE(pr.sources[1].kind == SourceKind::RainfallRate);
  REQUIRE(pr.sources[1].rate == Approx(1.0e-5).margin(1e-15));
  REQUIRE(pr.sources[2].kind == SourceKind::ExtractionWell);
  REQUIRE(pr.sources[2].rate == Approx(0.05).margin(1e-12));

  // The polygon contains check is wired through (sanity on a parsed ring).
  REQUIRE(pr.boundaries[0].polygon.contains(9.5, 5.0));
  REQUIRE_FALSE(pr.boundaries[0].polygon.contains(1.0, 5.0));
}

TEST_CASE("yaml polygon: absent sections => empty (b1-sw / b2-gw safe)") {
  const char* bare = R"(schema_version: '2.0'
simulation: {id: bare, mode: surface_water}
domain: {nx: 4, ny: 4, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: true, groundwater: false, solute: false}
sources: {surface: {rainfall: {from_file: false, rate: 0.0}}}
)";
  Config cfg = Config::fromString(bare, "");
  PolygonRegions pr = parsePolygonRegions(cfg);
  REQUIRE(pr.boundaries.empty());
  REQUIRE(pr.sources.empty());  // scalar `sources:` map is NOT a polygon sequence
}

TEST_CASE("yaml polygon: malformed entries fail loud") {
  const char* bad_type = R"(schema_version: '2.0'
simulation: {id: x, mode: surface_water}
domain: {nx: 4, ny: 4, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: true, groundwater: false, solute: false}
boundaries:
  - {name: bogus, type: bc_nonsense, vertices: [[0,0],[1,0],[1,1],[0,1]]}
)";
  bool threw = false;
  try {
    parsePolygonRegions(Config::fromString(bad_type, ""));
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);

  const char* few_verts = R"(schema_version: '2.0'
simulation: {id: x, mode: surface_water}
domain: {nx: 4, ny: 4, nz: 1, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: true, groundwater: false, solute: false}
sources:
  - {name: line, type: inflow_rate, vertices: [[0,0],[1,1]], rate: 1.0}
)";
  threw = false;
  try {
    parsePolygonRegions(Config::fromString(few_verts, ""));
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);
}
