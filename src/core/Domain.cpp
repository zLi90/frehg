#include "core/Domain.hpp"

namespace frehg2 {

Domain::Domain(const Grid& grid)
    : z("Domain::z", grid.nSurfaceCell()),
      area("Domain::area", grid.nSurfaceCell()),
      actMask("Domain::actMask", grid.nSurfaceCell()),
      roughness("Domain::roughness", grid.nSurfaceCell()),
      grid_(grid) {
  // Default: all cells active, unit area, zero bed, zero roughness (overwritten by IC/IO).
  Kokkos::deep_copy(actMask, 1.0);
  Kokkos::deep_copy(area, grid.dx() * grid.dy());
  Kokkos::deep_copy(z, 0.0);
  Kokkos::deep_copy(roughness, 0.0);
}

}  // namespace frehg2
