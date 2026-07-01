// P4.1 acceptance: SWE initialization (h = max(0, eta - z), u = v = 0) for flat and sloped
// beds, plus the Config -> SweSolver path on the frozen benchmark YAMLs (P4.0.1).
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include <string>

#include "core/Grid.hpp"
#include "frehg2_test.hpp"
#include "io/Config.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

namespace {
const char* kBench = FREHG2_BENCH_DIR;
std::string bench(const std::string& rel) { return std::string(kBench) + "/" + rel; }
}  // namespace

TEST_CASE("flat bed z=0, init_eta=2.0 -> h=2.0 everywhere, u=v=0") {
  Grid grid(5, 5, 1, 1.0, 1.0, 1.0);
  SweSolver swe(grid);
  SweParams p;
  p.min_depth = 1.0e-8;
  swe.setParams(p);
  swe.setBathymetryConstant(0.0);
  swe.initializeState(2.0);

  const auto& f = swe.fields();
  for (int j = 0; j < 5; ++j) {
    for (int i = 0; i < 5; ++i) {
      const int c = grid.getSurfaceIndex(i, j);
      REQUIRE(f.dept(c) == Approx(2.0).margin(1e-12));
      REQUIRE(f.eta(c) == Approx(2.0).margin(1e-12));
      REQUIRE(f.uu(c) == Approx(0.0).margin(1e-15));
      REQUIRE(f.vv(c) == Approx(0.0).margin(1e-15));
    }
  }
}

TEST_CASE("sloped bed z=j, init_eta=2.0 -> h=max(0,2-z), eta clamped to bed") {
  Grid grid(5, 5, 1, 1.0, 1.0, 1.0);
  SweSolver swe(grid);
  SweParams p;
  p.min_depth = 1.0e-8;
  swe.setParams(p);

  // Bathymetry rises with j: z(i,j) = j  (0..4).
  RealArr1DHost bed("bed", 25);
  for (int j = 0; j < 5; ++j)
    for (int i = 0; i < 5; ++i) bed(static_cast<size_t>(i + j * 5)) = static_cast<real>(j);
  swe.setBathymetry(bed);
  swe.initializeState(2.0);

  const double expect_dept[5] = {2.0, 1.0, 0.0, 0.0, 0.0};
  const auto& f = swe.fields();
  for (int j = 0; j < 5; ++j) {
    for (int i = 0; i < 5; ++i) {
      const int c = grid.getSurfaceIndex(i, j);
      REQUIRE(f.dept(c) == Approx(expect_dept[j]).margin(1e-12));
      // eta = max(init_eta, bed): clamped where bed > init_eta.
      const double expect_eta = (static_cast<double>(j) > 2.0) ? static_cast<double>(j) : 2.0;
      REQUIRE(f.eta(c) == Approx(expect_eta).margin(1e-12));
      REQUIRE(f.uu(c) == Approx(0.0).margin(1e-15));
      REQUIRE(f.vv(c) == Approx(0.0).margin(1e-15));
    }
  }
}

TEST_CASE("Config -> SweSolver path on b0-lake.yaml (parse + init)") {
  Config cfg = Config::fromFile(bench("b0-lake/b0-lake.yaml"));
  REQUIRE(cfg.get<int>("domain.nx") == 20);
  REQUIRE(cfg.get<bool>("modules.surface_water") == true);

  Grid grid(cfg.get<int>("domain.nx"), cfg.get<int>("domain.ny"), cfg.get<int>("domain.nz"),
            cfg.get<double>("domain.dx"), cfg.get<double>("domain.dy"),
            cfg.get<double>("domain.dz"));
  SweSolver swe(grid);
  swe.initialize(cfg);  // constant bed botz=0, init_eta=1.0
  REQUIRE(swe.params().gravity == Approx(9.81).margin(1e-12));
  REQUIRE(swe.params().manning == Approx(0.02).margin(1e-12));

  const auto& f = swe.fields();
  for (int j = 0; j < 20; ++j) {
    for (int i = 0; i < 20; ++i) {
      const int c = grid.getSurfaceIndex(i, j);
      REQUIRE(f.dept(c) == Approx(1.0).margin(1e-12));  // eta=1, bed=0
    }
  }
}

TEST_CASE("b1-sw.yaml parses and initializes (constant-bed smoke)") {
  Config cfg = Config::fromFile(bench("b1-sw/b1-sw.yaml"));
  Grid grid(cfg.get<int>("domain.nx"), cfg.get<int>("domain.ny"), cfg.get<int>("domain.nz"),
            cfg.get<double>("domain.dx"), cfg.get<double>("domain.dy"),
            cfg.get<double>("domain.dz"));
  SweSolver swe(grid);
  swe.initialize(cfg);
  REQUIRE(swe.params().min_depth == Approx(1.0e-8).margin(1e-20));
  // init_eta=-2.0, botz=-3.0 -> dept = 1.0 on the constant-bed smoke (real bath in b1 gate).
  const int c = grid.getSurfaceIndex(0, 4);
  REQUIRE(swe.fields().dept(c) == Approx(1.0).margin(1e-12));
}
