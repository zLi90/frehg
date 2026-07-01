// P8.3.1: solute concentration fields exist in the P2 state classes, are halo-sized, and
// are default-initialized to zero.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <array>
#include <vector>

#include <Kokkos_Core.hpp>

#include "core/Grid.hpp"
#include "core/GwState.hpp"
#include "core/State.hpp"

using namespace frehg2;

namespace {

bool allZero(const RealArr1D& v) {
  auto h = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h, v);
  for (size_t i = 0; i < h.extent(0); ++i)
    if (h(i) != 0.0) return false;
  return true;
}

}  // namespace

TEST_CASE("solute state: surface + subsurface conc exist, halo-sized, zero") {
  const std::vector<std::array<int, 3>> cases = {{4, 3, 5}, {7, 6, 2}};
  for (const auto& d : cases) {
    Grid grid(d[0], d[1], d[2], 1.0, 1.0, 0.1);
    State sw(grid);
    GwState gw(grid);

    // Surface conc is halo-padded (Grid::getSurfaceIndex addressable).
    REQUIRE(static_cast<int>(sw.conc.extent(0)) == grid.nSurfaceStorageCell());
    // Subsurface conc is halo-padded in all three dimensions.
    REQUIRE(static_cast<int>(gw.conc.extent(0)) == grid.nCell());

    REQUIRE(allZero(sw.conc));
    REQUIRE(allZero(gw.conc));
  }
}
