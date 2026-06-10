#include "swe/SweSolver.hpp"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <cmath>
#include <vector>

int main(int argc, char** argv)
{
#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, nullptr, nullptr);
#else
    (void)argc;
    (void)argv;
#endif
    frehg2::GridSpec spec;
    spec.nx = 4;
    spec.ny = 1;
    spec.nz = 1;

    frehg2::SweParameters params;
    params.dt = 10.0;
    params.rain_rate = 0.01;
    params.evap_rate = 0.002;
    params.min_depth = 1.0e-8;

    const frehg2::SweSolver solver(frehg2::Grid(spec), params);
    std::vector<frehg2::real> eta(4, 1.0);
    std::vector<frehg2::real> bed(4, 0.0);

    solver.applyRainEvap(eta, bed);
    for (const auto value : eta) {
        if (std::abs(value - 1.08) > 1.0e-12) {
#ifdef USE_PETSC
            PetscFinalize();
#endif
            return 1;
        }
    }

    params.rain_rate = 0.0;
    params.evap_rate = 1.0;
    const frehg2::SweSolver drying_solver(frehg2::Grid(spec), params);
    drying_solver.applyRainEvap(eta, bed);
    for (const auto value : eta) {
        if (value != 0.0) {
#ifdef USE_PETSC
            PetscFinalize();
#endif
            return 1;
        }
    }

    frehg2::GridSpec runoff_spec;
    runoff_spec.nx = 3;
    runoff_spec.ny = 1;
    runoff_spec.nz = 1;
    runoff_spec.dx = 1.0;
    runoff_spec.dy = 1.0;
    frehg2::SweParameters runoff_params;
    runoff_params.dt = 1.0;
    runoff_params.min_depth = 1.0e-8;
    runoff_params.wtfh = 1.0e-8;
    runoff_params.rain_rate = 1.0e-4;
    runoff_params.free_outflow_boundaries.push_back(
        {{frehg2::Grid(runoff_spec).getSurfaceIndex(2, 0)}, 1.0, 0.0});
    const frehg2::SweSolver runoff_solver(frehg2::Grid(runoff_spec), runoff_params);
    auto runoff_state = runoff_solver.initializeLegacyState({0.0, -0.01, -0.02}, 0.0);
    runoff_solver.advanceLegacyStep(runoff_state, runoff_params.rain_rate, 0.0);
    bool retained_rain = false;
    for (frehg2::index_t idx = 0; idx < runoff_solver.grid().nSurfaceCell(); ++idx) {
        retained_rain = retained_rain || runoff_state.dept[idx] > 0.0;
    }
    if (!retained_rain) {
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
