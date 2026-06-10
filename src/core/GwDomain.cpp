#include "core/GwDomain.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <utility>

namespace frehg2 {

GwDomain::GwDomain(Grid grid)
    : grid_(std::move(grid))
#ifdef USE_KOKKOS
      ,
      soilID("groundwater_soil_id", grid_.nCellMem()),
      dz("groundwater_layer_thickness", grid_.nz())
#endif
{
#ifdef USE_KOKKOS
    auto dz_h = Kokkos::create_mirror_view(dz);
    for (int k = 0; k < grid_.nz(); ++k) {
        dz_h(static_cast<std::size_t>(k)) = grid_.dz(k);
    }
    Kokkos::deep_copy(dz, dz_h);
    setUniformSoil(0);
#endif
}

const Grid& GwDomain::grid() const noexcept
{
    return grid_;
}

int GwDomain::nzGlob() const noexcept
{
    return grid_.nz();
}

real GwDomain::dzMultiplier() const noexcept
{
    return grid_.spec().dz_multiplier;
}

index_t GwDomain::size() const noexcept
{
    return grid_.nCellMem();
}

void GwDomain::setUniformSoil(int soil_id)
{
#ifdef USE_KOKKOS
    Kokkos::deep_copy(soilID, soil_id);
#else
    (void)soil_id;
#endif
}

}  // namespace frehg2
