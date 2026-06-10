#include "re/ReSolver.hpp"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <cmath>
#include <iostream>

int main()
{
#ifdef USE_PETSC
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);
#endif
    frehg2::GridSpec spec;
    spec.nx = 1;
    spec.ny = 1;
    spec.nz = 10;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.01;

    frehg2::ReParameters p;
    p.dt = 1.0e-4;
    p.dt_min = 1.0e-4;
    p.dt_max = 2.0;
    p.ksz = 2.89e-6;
    p.ss = 1.0e-5;
    p.soil_a = 1.43;
    p.soil_n = 1.56;
    p.wcs = 0.33;
    p.wcr = 0.0;
    p.init_wc = 0.033;
    p.htop = 0.0;
    p.hbot = 0.0;
    p.bc_type = {0, 0, 0, 0, 0, 1};
    p.dt_adjust = false;

    frehg2::ReSolver solver(frehg2::Grid(spec), p);
    auto state = solver.initializeLegacyState(0.0, -0.1);
    solver.computeConductivityFaces(state);
    auto predictor = solver.assemblePredictorSystem(state);
    const auto corrector_h = solver.solveLinearSystem(predictor);
    if (predictor.n != 10 || corrector_h.size() != 10) {
        std::cerr << "1D RE system size mismatch\n";
        return 1;
    }
    if (!(corrector_h.front() > state.h.front())) {
        std::cerr << "fixed top head should increase the top-cell pressure head\n";
        return 1;
    }

    state.h = corrector_h;
    solver.computeDarcyFlux(state);
    solver.updateWaterContent(state);
    if (!(state.wc.front() > p.init_wc)) {
        std::cerr << "top cell should wet under fixed head boundary\n";
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
