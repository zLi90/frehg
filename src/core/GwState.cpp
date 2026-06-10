#include "core/GwState.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <utility>

namespace frehg2 {

GwState::GwState(const Grid& grid)
    : size_(grid.nCellMem())
#ifdef USE_KOKKOS
      ,
      h_old("groundwater_head_old", size_),
      h_new("groundwater_head_new", size_),
      wc_old("groundwater_water_content_old", size_),
      wc_new("groundwater_water_content_new", size_),
      wc_excess("groundwater_water_content_excess", size_),
      k("groundwater_conductivity_faces", size_, 3),
      q("groundwater_darcy_flux_faces", size_, 3)
#endif
{
#ifdef USE_KOKKOS
    fill(0.0, 0.0);
    Kokkos::deep_copy(wc_excess, 0.0);
    Kokkos::deep_copy(k, 0.0);
    Kokkos::deep_copy(q, 0.0);
#endif
}

index_t GwState::size() const noexcept
{
    return size_;
}

void GwState::fill(real hydraulic_head, real water_content)
{
#ifdef USE_KOKKOS
    Kokkos::deep_copy(h_old, hydraulic_head);
    Kokkos::deep_copy(h_new, hydraulic_head);
    Kokkos::deep_copy(wc_old, water_content);
    Kokkos::deep_copy(wc_new, water_content);
#else
    (void)hydraulic_head;
    (void)water_content;
#endif
}

void GwState::swapTimeLevels()
{
#ifdef USE_KOKKOS
    std::swap(h_old, h_new);
    std::swap(wc_old, wc_new);
#endif
}

}  // namespace frehg2
