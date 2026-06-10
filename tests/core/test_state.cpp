#include "core/GwState.hpp"
#include "core/State.hpp"

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
        spec.nx = 3;
        spec.ny = 2;
        spec.nz = 4;
        const frehg2::Grid grid(spec);

        frehg2::State state(grid);
        frehg2::GwState gw_state(grid);

        if (state.size() != grid.nSurfaceCellMem()) {
            result = 1;
        }
        if (gw_state.size() != grid.nCellMem()) {
            result = 1;
        }

#ifdef USE_KOKKOS
        state.fill(1.5, 2.5, 3.5);
        Kokkos::deep_copy(state.h_new, 9.0);
        state.swapTimeLevels();

        auto h_old_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.h_old);
        auto h_new_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.h_new);
        auto hu_old_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.hu_old);

        if (h_old_h(0) != 9.0 || h_new_h(0) != 1.5 || hu_old_h(0) != 2.5) {
            result = 1;
        }
        if (state.qss.extent(0) != grid.nSurfaceCellMem()) {
            result = 1;
        }

        gw_state.fill(4.0, 0.35);
        Kokkos::deep_copy(gw_state.h_new, 8.0);
        Kokkos::deep_copy(gw_state.wc_new, 0.45);
        gw_state.swapTimeLevels();

        auto gw_h_old_h =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gw_state.h_old);
        auto gw_wc_old_h =
            Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gw_state.wc_old);

        if (gw_h_old_h(0) != 8.0 || gw_wc_old_h(0) != 0.45) {
            result = 1;
        }
        if (gw_state.k.extent(0) != grid.nCellMem() || gw_state.k.extent(1) != 3) {
            result = 1;
        }
        if (gw_state.q.extent(0) != grid.nCellMem() || gw_state.q.extent(1) != 3) {
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
