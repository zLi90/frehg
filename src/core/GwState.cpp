#include "core/GwState.hpp"

namespace frehg2 {

GwState::GwState(const Grid& grid)
    : h("GwState::h", grid.nActiveCell()),
      hn("GwState::hn", grid.nActiveCell()),
      wc("GwState::wc", grid.nActiveCell()),
      wcn("GwState::wcn", grid.nActiveCell()),
      Kx("GwState::Kx", grid.nActiveCell()),
      Ky("GwState::Ky", grid.nActiveCell()),
      Kz("GwState::Kz", grid.nActiveCell()),
      qx("GwState::qx", grid.nActiveCell()),
      qy("GwState::qy", grid.nActiveCell()),
      qz("GwState::qz", grid.nActiveCell()),
      h_pred("GwState::h_pred", grid.nActiveCell()),
      wc_pred("GwState::wc_pred", grid.nActiveCell()),
      conc("GwState::conc", grid.nCell()),
      grid_(grid) {
  Kokkos::deep_copy(h, 0.0);
  Kokkos::deep_copy(hn, 0.0);
  Kokkos::deep_copy(wc, 0.0);
  Kokkos::deep_copy(wcn, 0.0);
  Kokkos::deep_copy(Kx, 0.0);
  Kokkos::deep_copy(Ky, 0.0);
  Kokkos::deep_copy(Kz, 0.0);
  Kokkos::deep_copy(qx, 0.0);
  Kokkos::deep_copy(qy, 0.0);
  Kokkos::deep_copy(qz, 0.0);
  Kokkos::deep_copy(h_pred, 0.0);
  Kokkos::deep_copy(wc_pred, 0.0);
  Kokkos::deep_copy(conc, 0.0);
}

}  // namespace frehg2
