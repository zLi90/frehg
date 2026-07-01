// P12.3.2 acceptance: PolygonIndex on a 100x100 grid with one square polygon. Cells whose
// centroids fall inside are tagged; boundary cells just inside/outside are correct; the global
// count matches the analytic cell count.
#define FREHG2_TEST_IMPL
#include "bc/PolygonIndex.hpp"
#include "bc/Polygon.hpp"
#include "core/Grid.hpp"
#include "frehg2_test.hpp"

using namespace frehg2;

TEST_CASE("polygon index: square region on a 100x100 unit grid") {
  // Unit cells, origin (0,0): cell (i,j) centroid at (i+0.5, j+0.5).
  Grid g(100, 100, 1, 1.0, 1.0, 1.0);
  Polygon sq;
  sq.name = "block";
  // x,y in [20,40] => centroids 20.5..39.5 => i,j in [20,39] (20 cells each way).
  sq.vertices = {{20.0, 20.0}, {40.0, 20.0}, {40.0, 40.0}, {20.0, 40.0}};

  PolygonIndex idx;
  idx.build({sq}, g, nullptr, 0.0, 0.0);
  REQUIRE(idx.nPolygons() == 1);

  // Interior cell tagged.
  REQUIRE(idx.lookup(g.getSurfaceIndex(30, 30)) == 0);
  // Just-inside boundary cells.
  REQUIRE(idx.lookup(g.getSurfaceIndex(20, 20)) == 0);
  REQUIRE(idx.lookup(g.getSurfaceIndex(39, 39)) == 0);
  // Just-outside cells.
  REQUIRE(idx.lookup(g.getSurfaceIndex(19, 30)) == -1);
  REQUIRE(idx.lookup(g.getSurfaceIndex(40, 30)) == -1);
  REQUIRE(idx.lookup(g.getSurfaceIndex(30, 19)) == -1);
  REQUIRE(idx.lookup(g.getSurfaceIndex(30, 40)) == -1);
  // A far cell.
  REQUIRE(idx.lookup(g.getSurfaceIndex(90, 90)) == -1);

  // Analytic count: 20 x 20 = 400 cells.
  REQUIRE(idx.countInPolygon(0) == 400);
  auto counts = idx.globalColumnCounts(nullptr);
  REQUIRE(counts.size() == 1u);
  REQUIRE(counts[0] == 400);
}

TEST_CASE("polygon index: first-match precedence on overlap; empty list tags nothing") {
  Grid g(10, 10, 1, 1.0, 1.0, 1.0);
  Polygon a;
  a.name = "a";
  a.vertices = {{0.0, 0.0}, {6.0, 0.0}, {6.0, 6.0}, {0.0, 6.0}};
  Polygon b;
  b.name = "b";
  b.vertices = {{4.0, 4.0}, {10.0, 4.0}, {10.0, 10.0}, {4.0, 10.0}};

  PolygonIndex idx;
  idx.build({a, b}, g, nullptr, 0.0, 0.0);
  // Cell (5,5) is in both a and b; first (index 0) wins.
  REQUIRE(idx.lookup(g.getSurfaceIndex(5, 5)) == 0);
  REQUIRE(idx.lookup(g.getSurfaceIndex(1, 1)) == 0);
  REQUIRE(idx.lookup(g.getSurfaceIndex(8, 8)) == 1);

  PolygonIndex empty;
  empty.build({}, g, nullptr, 0.0, 0.0);
  REQUIRE(empty.empty());
  REQUIRE(empty.lookup(g.getSurfaceIndex(5, 5)) == -1);
}

TEST_CASE("polygon index: nonzero origin shifts the mapping") {
  Grid g(10, 10, 1, 2.0, 2.0, 1.0);  // dx=dy=2
  // With origin (100,200), cell (i,j) centroid = (100 + (i+0.5)*2, 200 + (j+0.5)*2).
  // A polygon around centroid of cell (3,4): (100+7, 200+9) = (107, 209).
  Polygon sq;
  sq.vertices = {{106.0, 208.0}, {108.0, 208.0}, {108.0, 210.0}, {106.0, 210.0}};
  PolygonIndex idx;
  idx.build({sq}, g, nullptr, 100.0, 200.0);
  REQUIRE(idx.lookup(g.getSurfaceIndex(3, 4)) == 0);
  REQUIRE(idx.countInPolygon(0) == 1);
}
