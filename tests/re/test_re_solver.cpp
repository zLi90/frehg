#include "re/ReSolver.hpp"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <iostream>

int main()
{
#ifdef USE_PETSC
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
#endif
    frehg2::GridSpec spec;
    spec.nx = 1;
    spec.ny = 1;
    spec.nz = 20;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.01;

    frehg2::ReParameters p;
    p.dt = 1.0e-4;
    p.ksz = 2.89e-6;
    p.ss = 1.0e-5;
    p.soil_a = 1.43;
    p.soil_n = 1.56;
    p.wcs = 0.33;
    p.wcr = 0.0;
    p.init_wc = 0.033;
    p.htop = 0.0;
    p.bc_type = {0, 0, 0, 0, 0, 1};
    p.dt_adjust = true;

    frehg2::ReSolver solver(frehg2::Grid(spec), p);
    auto state = solver.initializeLegacyState(0.0, -0.2);
    const auto h0 = state.h.front();
    solver.advanceLegacyStep(state);
    if (!(state.h.front() > h0)) {
        std::cerr << "RE step did not advance hydraulic head\n";
#ifdef USE_PETSC
        PetscFinalize();
#endif
        return 1;
    }
    if (state.dtg < p.dt_min || state.dtg > p.dt_max) {
        std::cerr << "adaptive groundwater timestep outside bounds\n";
#ifdef USE_PETSC
        PetscFinalize();
#endif
        return 1;
    }
#ifdef USE_PETSC
    PetscFinalize();
#endif
    return 0;
}
