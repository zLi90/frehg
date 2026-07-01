// P12.3.1 acceptance: Polygon::contains() ray-casting on square, concave (L-shape), and a
// donut (keyhole ring with an interior hole). Known inside/outside points are asserted; the
// bounding-box fast-reject is exercised by points outside the extent.
#define FREHG2_TEST_IMPL
#include "bc/Polygon.hpp"
#include "frehg2_test.hpp"

using namespace frehg2;

TEST_CASE("polygon: axis-aligned square") {
  Polygon sq;
  sq.name = "square";
  sq.vertices = {{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};

  REQUIRE(sq.contains(5.0, 5.0));
  REQUIRE(sq.contains(0.01, 0.01));
  REQUIRE(sq.contains(9.99, 9.99));
  REQUIRE_FALSE(sq.contains(-0.5, 5.0));
  REQUIRE_FALSE(sq.contains(10.5, 5.0));
  REQUIRE_FALSE(sq.contains(5.0, -0.5));
  REQUIRE_FALSE(sq.contains(5.0, 10.5));
  // Far points (bounding-box reject).
  REQUIRE_FALSE(sq.contains(1000.0, 1000.0));
}

TEST_CASE("polygon: concave L-shape") {
  // L-shape occupying the bottom row (y<4) full width and the left column (x<4) full height.
  Polygon ell;
  ell.name = "L";
  ell.vertices = {{0.0, 0.0}, {10.0, 0.0}, {10.0, 4.0}, {4.0, 4.0},
                  {4.0, 10.0}, {0.0, 10.0}};

  REQUIRE(ell.contains(2.0, 2.0));   // lower-left corner: in
  REQUIRE(ell.contains(8.0, 2.0));   // lower-right arm: in
  REQUIRE(ell.contains(2.0, 8.0));   // upper-left arm: in
  REQUIRE_FALSE(ell.contains(8.0, 8.0));  // notch (concave cut): out
  REQUIRE_FALSE(ell.contains(5.0, 5.0));  // notch: out
  REQUIRE(ell.contains(3.9, 9.0));   // just inside the left arm
  REQUIRE_FALSE(ell.contains(4.1, 9.0));  // just outside the left arm (in the notch)
}

TEST_CASE("polygon: donut (keyhole ring with interior hole)") {
  // Outer 10x10 square with a 4..6 square hole, expressed as a single self-touching keyhole ring
  // (seam on the left edge at y=5). Even-odd rule makes the hole interior "outside".
  Polygon donut;
  donut.name = "donut";
  donut.vertices = {{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}, {0.0, 5.0},
                    {4.0, 5.0}, {4.0, 4.0},  {6.0, 4.0},   {6.0, 6.0},  {4.0, 6.0},
                    {4.0, 5.0}, {0.0, 5.0}};

  // Material (the ring) is inside.
  REQUIRE(donut.contains(1.0, 5.0));
  REQUIRE(donut.contains(8.0, 5.0));
  REQUIRE(donut.contains(5.0, 8.0));
  REQUIRE(donut.contains(5.0, 2.0));
  REQUIRE(donut.contains(1.0, 1.0));
  REQUIRE(donut.contains(9.0, 9.0));
  // The hole is outside.
  REQUIRE_FALSE(donut.contains(5.0, 5.0));
  REQUIRE_FALSE(donut.contains(4.5, 5.0));
  REQUIRE_FALSE(donut.contains(5.0, 5.5));
  // Outside the outer square.
  REQUIRE_FALSE(donut.contains(-1.0, 5.0));
  REQUIRE_FALSE(donut.contains(11.0, 5.0));
}

TEST_CASE("polygon: degenerate rings never contain") {
  Polygon empty;
  REQUIRE_FALSE(empty.contains(0.0, 0.0));
  Polygon seg;
  seg.vertices = {{0.0, 0.0}, {1.0, 1.0}};
  REQUIRE_FALSE(seg.contains(0.5, 0.5));
}
