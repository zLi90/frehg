#include "core/Domain.hpp"
#include "core/State.hpp"
#include "swe/SweSolver.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

int main(int argc, char** argv)
{
#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    int result = 0;
    try {
        frehg2::GridSpec spec;
        spec.nx = 3;
        spec.ny = 2;
        spec.nz = 1;
        frehg2::Grid grid(spec);
        frehg2::Domain domain(grid);
        frehg2::State state(grid);
        domain.setUniformBed(2.0);

        frehg2::SweSolver solver(grid, frehg2::SweParameters{});
        solver.initialize(domain, state, 5.5);

#ifdef USE_KOKKOS
        const auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.h_old);
        const auto z = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), state.z);
        if (h(0) != 3.5 || z(0) != 2.0) {
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
