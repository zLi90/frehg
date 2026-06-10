#include "core/Domain.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <utility>

namespace frehg2 {

Domain::Domain(Grid grid)
    : grid_(std::move(grid))
#ifdef USE_KOKKOS
      ,
      z("surface_bed_elevation", grid_.nSurfaceCellMem()),
      area("surface_cell_area", grid_.nSurfaceCellMem()),
      actMask("surface_active_mask", grid_.nSurfaceCellMem())
#endif
{
#ifdef USE_KOKKOS
    Kokkos::deep_copy(area, grid_.dx() * grid_.dy());
    setAllActive(1);
#endif
}

const Grid& Domain::grid() const noexcept
{
    return grid_;
}

index_t Domain::size() const noexcept
{
    return grid_.nSurfaceCellMem();
}

void Domain::setUniformBed(real bed_elevation)
{
#ifdef USE_KOKKOS
    Kokkos::deep_copy(z, bed_elevation);
#else
    (void)bed_elevation;
#endif
}

void Domain::setAllActive(int active)
{
#ifdef USE_KOKKOS
    Kokkos::deep_copy(actMask, active);
#else
    (void)active;
#endif
}

}  // namespace frehg2
