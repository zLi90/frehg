// P2.1 acceptance: Grid halo-padded indexing + LegacyIndexAdapter round-trip.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include "core/Grid.hpp"
#include "frehg2/core/LegacyIndexAdapter.hpp"

using namespace frehg2;

TEST_CASE("surface/cell counts") {
  Grid g(10, 10, 1, 1.0, 1.0, 1.0);
  REQUIRE(g.nSurfaceCell() == 100);

  Grid g3(10, 10, 5, 1.0, 1.0, 1.0);
  REQUIRE(g3.nCell() == 12 * 12 * 7);   // horizontal + vertical halo
  REQUIRE(g3.nCell() == 1008);
  REQUIRE(g3.nActiveCell() == 10 * 10 * 5);
  REQUIRE(g3.nSurfaceStorageCell() == 12 * 12);
}

TEST_CASE("getIndex matches the halo formula") {
  Grid g(10, 10, 5, 1.0, 1.0, 1.0);
  const int nx = 10;
  REQUIRE(g.getIndex(0, 0, 0) ==
          (0 + 1) + (0 + 1) * (nx + 2) + (0 + 1) * (nx + 2) * (nx + 2));
  REQUIRE(g.getIndex(0, 0, 0) == 157);
  REQUIRE(g.getSurfaceIndex(0, 0) == 13);
  REQUIRE(g.getIndex(0, 0, -1) == 13);
}

TEST_CASE("getIndex/getIJK round-trip for all storage cells in 4x4x3") {
  Grid g(4, 4, 3, 1.0, 1.0, 1.0);
  for (int flat = 0; flat < g.nCell(); ++flat) {
    int i = 0, j = 0, k = 0;
    g.getIJK(flat, i, j, k);
    REQUIRE(g.getIndex(i, j, k) == flat);
  }
}

TEST_CASE("halo cells are inactive") {
  Grid g(4, 4, 3, 1.0, 1.0, 1.0);
  REQUIRE_FALSE(g.isActive(-1, 0, 0));
  REQUIRE_FALSE(g.isActive(4, 0, 0));
  REQUIRE_FALSE(g.isActive(0, -1, 0));
  REQUIRE_FALSE(g.isActive(0, 4, 0));
  REQUIRE(g.isActive(0, 0, 0));
  REQUIRE(g.isActive(3, 3, 2));
  REQUIRE_FALSE(g.isActive(0, 0, -1));
  REQUIRE_FALSE(g.isActive(0, 0, 3));
}

TEST_CASE("LegacyIndexAdapter index mapping") {
  // legacy index 5 in a 10x10 grid -> (5%10, 5/10) = (5, 0)
  REQUIRE(LegacyIndexAdapter::legacySurfaceInteriorIndex(5, 0, 10) == 5);
  // Frehg2 coordinate (5,0) -> surface storage (5+1) + (0+1)*(10+2) = 18
  REQUIRE(LegacyIndexAdapter::frehg2SurfaceStorageIndex(5, 0, 10) == 18);
  // subsurface: (i+j*nx)*nz + k
  REQUIRE(LegacyIndexAdapter::legacySubsurfaceInteriorIndex(2, 1, 3, 4, 5) ==
          (2 + 1 * 4) * 5 + 3);
}

TEST_CASE("LegacyIndexAdapter surface order round-trip") {
  const int nx = 6, ny = 4;
  Grid g(nx, ny, 1, 1.0, 1.0, 1.0);
  RealArr1DHost legacy("legacy", static_cast<size_t>(nx) * ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) legacy(i + j * nx) = static_cast<real>(100 * j + i);

  auto halo = LegacyIndexAdapter::fromLegacySurfaceOrder(legacy, nx, ny);
  // ghost cells must be zero
  REQUIRE(halo(g.getSurfaceIndex(-1, 0)) == Approx(0.0).margin(1e-15));
  // interior maps to the halo formula
  REQUIRE(halo(g.getSurfaceIndex(3, 2)) == Approx(100 * 2 + 3).margin(1e-12));

  auto back = LegacyIndexAdapter::toLegacySurfaceOrder(halo, nx, ny);
  for (int idx = 0; idx < nx * ny; ++idx)
    REQUIRE(back(idx) == Approx(legacy(idx)).margin(1e-12));
}

TEST_CASE("LegacyIndexAdapter subsurface order round-trip") {
  const int nx = 3, ny = 2, nz = 4;
  RealArr1DHost legacy("legacy3d", static_cast<size_t>(nx) * ny * nz);
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        legacy((i + j * nx) * nz + k) = static_cast<real>(1000 * k + 10 * j + i);

  auto halo = LegacyIndexAdapter::fromLegacySubsurfaceOrder(legacy, nx, ny, nz);
  auto back = LegacyIndexAdapter::toLegacySubsurfaceOrder(halo, nx, ny, nz);
  for (int idx = 0; idx < nx * ny * nz; ++idx)
    REQUIRE(back(idx) == Approx(legacy(idx)).margin(1e-12));
}
