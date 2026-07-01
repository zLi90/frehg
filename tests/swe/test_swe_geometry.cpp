// P4 stage 4a: SWE face geometry (Asx/Asy/Asz/Vs/Vsx/Vsy/Aszx/Aszy) on hand-computable
// cases. This is the unit-level geometry check; element-wise parity against the
// instrumented legacy dump (<1e-12) is the subsequent 4b gate (test_swe_geometry_legacy).
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "core/Grid.hpp"
#include "frehg2_test.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

TEST_CASE("flat wet bed geometry: Asz=dx*dy, Asx=deptx*dy, Asy=depty*dx, Vs=dept*dx*dy") {
  const real dx = 2.0, dy = 3.0;
  Grid grid(3, 3, 1, dx, dy, 1.0);
  SweSolver swe(grid);
  SweParams p;
  p.min_depth = 1.0e-9;
  swe.setParams(p);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(1.0);  // eta = 1 -> dept = 1 everywhere

  const auto& f = swe.fields();
  const int c = grid.getSurfaceIndex(1, 1);  // center
  REQUIRE(f.dept(c) == Approx(1.0).margin(1e-12));
  REQUIRE(f.Asz(c) == Approx(dx * dy).margin(1e-12));        // 6
  REQUIRE(f.Asx(c) == Approx(1.0 * dy).margin(1e-12));       // deptx*dy = 3
  REQUIRE(f.Asy(c) == Approx(1.0 * dx).margin(1e-12));       // depty*dx = 2
  REQUIRE(f.Vs(c) == Approx(1.0 * dx * dy).margin(1e-12));   // 6
  REQUIRE(f.Vsx(c) == Approx(dx * dy).margin(1e-12));        // 0.5*(6+6) = 6
  REQUIRE(f.Vsy(c) == Approx(dx * dy).margin(1e-12));
  REQUIRE(f.Aszx(c) == Approx(dx * dy).margin(1e-12));
  REQUIRE(f.Aszy(c) == Approx(dx * dy).margin(1e-12));
}

TEST_CASE("dry cell geometry: Asz=0, face areas vanish across the dry face") {
  const real dx = 1.0, dy = 1.0;
  Grid grid(3, 3, 1, dx, dy, 1.0);
  SweSolver swe(grid);
  SweParams p;
  p.min_depth = 1.0e-9;
  swe.setParams(p);

  // Center cell bed is high and stays dry; the rest is wet at eta=1.
  RealArr1DHost bed("bed", 9);
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i) bed(static_cast<size_t>(i + j * 3)) = 0.0;
  bed(static_cast<size_t>(1 + 1 * 3)) = 5.0;  // center
  swe.setBathymetry(bed);
  swe.initializeState(1.0);

  const auto& f = swe.fields();
  const int center = grid.getSurfaceIndex(1, 1);
  const int left = grid.getSurfaceIndex(0, 1);
  REQUIRE(f.dept(center) == Approx(0.0).margin(1e-12));   // 5-5 clamp -> 0
  REQUIRE(f.Asz(center) == Approx(0.0).margin(1e-12));
  REQUIRE(f.dept(left) == Approx(1.0).margin(1e-12));
  REQUIRE(f.Asz(left) == Approx(dx * dy).margin(1e-12));
  // x+ face of the left cell borders the dry/high center: deptx = max(eta)-max(bed)
  // = max(1,5) - max(0,5) = 0 -> Asx = 0.
  REQUIRE(f.deptx(left) == Approx(0.0).margin(1e-12));
  REQUIRE(f.Asx(left) == Approx(0.0).margin(1e-12));
  // Aszx across the wet/dry face = 0.5*(Asz_left + Asz_center) = 0.5*(1 + 0) = 0.5.
  REQUIRE(f.Aszx(left) == Approx(0.5 * dx * dy).margin(1e-12));
}

TEST_CASE("CFL diagnostic formula: still water has cfl = sqrt(g*h)*dt/min(dx,dy)") {
  Grid grid(4, 4, 1, 10.0, 10.0, 1.0);
  SweSolver swe(grid);
  SweParams p;
  p.gravity = 9.81;
  p.dt = 2.0;
  p.min_depth = 1.0e-9;
  swe.setParams(p);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(4.0);  // h = 4, u = v = 0

  const double expect = std::sqrt(9.81 * 4.0) * 2.0 / 10.0;
  REQUIRE(swe.maxCfl() == Approx(expect).margin(1e-12));
}
