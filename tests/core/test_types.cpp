// P1.3 acceptance: core type definitions compile and have the expected properties.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <type_traits>

#include "frehg2/core/define.hpp"
#include "frehg2/core/types.hpp"

using namespace frehg2;

// Compile-time guarantees.
static_assert(sizeof(real) == 8, "real must be double precision (8 bytes)");
static_assert(std::is_same<real, double>::value, "real must be double");
static_assert(sizeof(index_t) == 4, "index_t must be 32-bit");

TEST_CASE("real is double precision") {
  REQUIRE(sizeof(real) == 8);
  REQUIRE(static_cast<real>(0.1) == Approx(0.1).margin(1e-15));
}

TEST_CASE("NODATA sentinel and enums") {
  REQUIRE(NO_DATA_REAL == Approx(-9999.0).margin(1e-12));
  REQUIRE(static_cast<int>(SimMode::SW_ONLY) == 0);
  REQUIRE(static_cast<int>(SimMode::COUPLED) == 2);
  REQUIRE(static_cast<int>(BCType::FIXED_HEAD) == 4);
}

TEST_CASE("device and host Kokkos views allocate with correct extents") {
  RealArr1D a("a", 10);
  REQUIRE(a.extent(0) == 10u);

  RealArr3D g("g", 2, 3, 4);
  REQUIRE(g.extent(0) == 2u);
  REQUIRE(g.extent(1) == 3u);
  REQUIRE(g.extent(2) == 4u);

  IntArr1DHost h("h", 5);
  REQUIRE(h.extent(0) == 5u);

  // Host mirror round-trip of a device view.
  RealArr1D dev("dev", 4);
  auto mirror = Kokkos::create_mirror_view(dev);
  for (int i = 0; i < 4; ++i) {
    mirror(i) = static_cast<real>(i) * 2.0;
  }
  Kokkos::deep_copy(dev, mirror);
  Kokkos::deep_copy(mirror, dev);
  REQUIRE(mirror(3) == Approx(6.0).margin(1e-12));
}

TEST_CASE("parameter structs hold values with two distinct cases") {
  DomainParams d1{4, 5, 6, 1.0, 2.0, 0.1, 1.05, -3.0};
  REQUIRE(d1.nx == 4);
  REQUIRE(d1.nz == 6);
  REQUIRE(d1.dz_incre == Approx(1.05).margin(1e-12));

  // Second parameter set (edge-ish: single column, geometric off).
  DomainParams d2{1, 1, 100, 0.01, 0.01, 0.01, 1.0, -1.0};
  REQUIRE(d2.nz == 100);
  REQUIRE(d2.bot_z == Approx(-1.0).margin(1e-12));

  TimeParams t{5.0, 18000.0, 1800.0, 100, 0.0, 0};
  REQUIRE(t.t_end == Approx(18000.0).margin(1e-9));
  REQUIRE(t.max_steps == 100);

  ModuleFlags m{true, false, false};
  REQUIRE(m.surface_water);
  REQUIRE_FALSE(m.groundwater);
}
