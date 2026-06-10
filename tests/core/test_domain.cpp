#include "core/Domain.hpp"
#include "core/GwDomain.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

int main(int argc, char** argv)
{
#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#endif

    int result = 0;

    try {
        frehg2::GridSpec spec;
        spec.nx = 4;
        spec.ny = 3;
        spec.nz = 2;
        spec.dx = 5.0;
        spec.dy = 6.0;
        spec.dz = 0.25;
        spec.dz_multiplier = 4.0;

        const frehg2::Grid grid(spec);
        frehg2::Domain domain(grid);
        frehg2::GwDomain gw_domain(grid);

        if (domain.size() != grid.nSurfaceCellMem()) {
            result = 1;
        }
        if (gw_domain.size() != grid.nCellMem()) {
            result = 1;
        }
        if (gw_domain.nzGlob() != 2 || gw_domain.dzMultiplier() != 4.0) {
            result = 1;
        }

#ifdef USE_KOKKOS
        if (domain.z.extent(0) != grid.nSurfaceCellMem()) {
            result = 1;
        }
        if (domain.area.extent(0) != grid.nSurfaceCellMem()) {
            result = 1;
        }
        if (domain.actMask.extent(0) != grid.nSurfaceCellMem()) {
            result = 1;
        }
        if (gw_domain.soilID.extent(0) != grid.nCellMem()) {
            result = 1;
        }
        if (gw_domain.dz.extent(0) != 2) {
            result = 1;
        }

        domain.setUniformBed(7.0);
        domain.setAllActive(1);
        gw_domain.setUniformSoil(3);

        const auto z_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), domain.z);
        const auto area_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), domain.area);
        const auto mask_h =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), domain.actMask);
        const auto soil_h =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gw_domain.soilID);
        const auto dz_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gw_domain.dz);

        if (z_h(0) != 7.0 || area_h(0) != 30.0 || mask_h(0) != 1 || soil_h(0) != 3) {
            result = 1;
        }
        if (dz_h(0) != 0.25 || dz_h(1) != 1.0) {
            result = 1;
        }
#endif
    } catch (...) {
        result = 1;
    }

#ifdef USE_KOKKOS
    Kokkos::finalize();
#endif

    return result;
}
