#include "core/GwDomain.hpp"

namespace frehg2 {

GwDomain::GwDomain(const Grid& grid, real bot_z)
    : dz3d("GwDomain::dz3d", grid.nz()),
      soilID("GwDomain::soilID", grid.nActiveCell()),
      z3d("GwDomain::z3d", grid.nActiveCell()),
      grid_(grid),
      bot_z_(bot_z) {
  const int nz = grid.nz();
  const real dz0 = grid.dz();
  const real ratio = grid.dzIncre();

  // dz3d[k] = dz * dz_incre^k (host build, then copy to device).
  auto dz_h = Kokkos::create_mirror_view(dz3d);
  real factor = 1.0;
  real total = 0.0;
  for (int k = 0; k < nz; ++k) {
    dz_h(k) = dz0 * factor;
    total += dz_h(k);
    factor *= ratio;
  }
  Kokkos::deep_copy(dz3d, dz_h);

  // z3d: absolute cell-center elevation, top (k=0) to bottom (k=nz-1). The domain top is
  // bot_z + total thickness; each column shares the same vertical profile here (refined by
  // per-column surface elevation in P5).
  const real top = bot_z + total;
  const int nx = grid.nx();
  const int ny = grid.ny();
  auto z3d_h = Kokkos::create_mirror_view(z3d);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      real depth_to_top_of_layer = 0.0;
      for (int k = 0; k < nz; ++k) {
        const real center = top - depth_to_top_of_layer - 0.5 * dz_h(k);
        z3d_h((i + j * nx) * nz + k) = center;
        depth_to_top_of_layer += dz_h(k);
      }
    }
  }
  Kokkos::deep_copy(z3d, z3d_h);

  Kokkos::deep_copy(soilID, 0);
}

}  // namespace frehg2
