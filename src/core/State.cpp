#include "core/State.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <utility>

namespace frehg2 {

State::State(const Grid& grid)
    : size_(grid.nSurfaceCellMem())
#ifdef USE_KOKKOS
      ,
      h_old("surface_water_depth_old", size_),
      h_new("surface_water_depth_new", size_),
      hu_old("surface_momentum_x_old", size_),
      hu_new("surface_momentum_x_new", size_),
      hv_old("surface_momentum_y_old", size_),
      hv_new("surface_momentum_y_new", size_),
      z("surface_bed_elevation_state", size_),
      roughness("surface_manning_roughness", size_),
      qss("surface_subsurface_exchange_flux", size_)
#endif
{
#ifdef USE_KOKKOS
    fill(0.0, 0.0, 0.0);
    Kokkos::deep_copy(roughness, 0.0);
    Kokkos::deep_copy(qss, 0.0);
    Kokkos::deep_copy(z, 0.0);
#endif
}

index_t State::size() const noexcept
{
    return size_;
}

void State::fill(real water_depth, real momentum_x, real momentum_y)
{
#ifdef USE_KOKKOS
    Kokkos::deep_copy(h_old, water_depth);
    Kokkos::deep_copy(h_new, water_depth);
    Kokkos::deep_copy(hu_old, momentum_x);
    Kokkos::deep_copy(hu_new, momentum_x);
    Kokkos::deep_copy(hv_old, momentum_y);
    Kokkos::deep_copy(hv_new, momentum_y);
#else
    (void)water_depth;
    (void)momentum_x;
    (void)momentum_y;
#endif
}

void State::swapTimeLevels()
{
#ifdef USE_KOKKOS
    std::swap(h_old, h_new);
    std::swap(hu_old, hu_new);
    std::swap(hv_old, hv_new);
#endif
}

}  // namespace frehg2
