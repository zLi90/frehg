// P2.2 acceptance: Domain / GwDomain.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include "core/Domain.hpp"
#include "core/GwDomain.hpp"

using namespace frehg2;

TEST_CASE("surface domain field sizes") {
  Grid g(10, 10, 1, 1.0, 1.0, 1.0);
  Domain d(g);
  REQUIRE(d.actMask.extent(0) == 100u);
  REQUIRE(d.z.extent(0) == 100u);
  REQUIRE(d.area.extent(0) == 100u);
  REQUIRE(d.roughness.extent(0) == 100u);

  auto mask = Kokkos::create_mirror_view(d.actMask);
  Kokkos::deep_copy(mask, d.actMask);
  REQUIRE(mask(0) == Approx(1.0).margin(1e-15));
}

TEST_CASE("geometric layer thicknesses dz3d = dz * dz_incre^k") {
  Grid g(2, 2, 5, 1.0, 1.0, 0.01, 1.1);
  GwDomain gd(g);
  REQUIRE(gd.dz3d.extent(0) == 5u);

  auto dz = Kokkos::create_mirror_view(gd.dz3d);
  Kokkos::deep_copy(dz, gd.dz3d);
  REQUIRE(dz(0) == Approx(0.01).margin(1e-12));
  REQUIRE(dz(4) == Approx(0.01 * 1.1 * 1.1 * 1.1 * 1.1).margin(1e-9));   // ~0.014641
  REQUIRE(dz(4) == Approx(0.0146410).margin(1e-7));

  // Second parameter set: uniform layers (dz_incre = 1).
  Grid g2(1, 1, 3, 1.0, 1.0, 0.5, 1.0);
  GwDomain gd2(g2);
  auto dz2 = Kokkos::create_mirror_view(gd2.dz3d);
  Kokkos::deep_copy(dz2, gd2.dz3d);
  REQUIRE(dz2(0) == Approx(0.5).margin(1e-12));
  REQUIRE(dz2(2) == Approx(0.5).margin(1e-12));
}

TEST_CASE("z3d is ordered top (k=0) to bottom (k=nz-1)") {
  Grid g(1, 1, 4, 1.0, 1.0, 0.25, 1.0);
  GwDomain gd(g, /*bot_z=*/-1.0);
  REQUIRE(gd.z3d.extent(0) == 4u);
  auto z = Kokkos::create_mirror_view(gd.z3d);
  Kokkos::deep_copy(z, gd.z3d);
  // top layer center must be higher than bottom layer center
  REQUIRE(z(0) > z(3));
  // total thickness = 4*0.25 = 1.0; top = bot_z + 1.0 = 0.0; layer0 center = -0.125
  REQUIRE(z(0) == Approx(-0.125).margin(1e-12));
  REQUIRE(z(3) == Approx(-0.875).margin(1e-12));
}

TEST_CASE("GW soilID/z3d sized to active cells") {
  Grid g(3, 2, 4, 1.0, 1.0, 0.1, 1.0);
  GwDomain gd(g);
  REQUIRE(gd.soilID.extent(0) == 24u);
  REQUIRE(gd.z3d.extent(0) == 24u);
}
