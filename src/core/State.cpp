#include "core/State.hpp"

namespace frehg2 {

State::State(const Grid& grid)
    : eta("State::eta", grid.nSurfaceCell()),
      u("State::u", grid.nSurfaceCell()),
      v("State::v", grid.nSurfaceCell()),
      h("State::h", grid.nSurfaceCell()),
      qss("State::qss", grid.nSurfaceCell()),
      conc("State::conc", grid.nSurfaceStorageCell()),
      grid_(grid) {
  Kokkos::deep_copy(eta, 0.0);
  Kokkos::deep_copy(u, 0.0);
  Kokkos::deep_copy(v, 0.0);
  Kokkos::deep_copy(h, 0.0);
  Kokkos::deep_copy(qss, 0.0);
  Kokkos::deep_copy(conc, 0.0);
}

}  // namespace frehg2
