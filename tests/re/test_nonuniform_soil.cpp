#include "re/ReSolver.hpp"

#include <cmath>
#include <iostream>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

}  // namespace

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 2;
    spec.ny = 1;
    spec.nz = 2;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.1;
    const frehg2::Grid grid(spec);

    frehg2::ReParameters p;
    p.init_wc = 0.20;
    p.soil_table = {
        frehg2::SoilParameters{0.0, 0.0, 1.0e-7, 1.0, 2.0, 0.40, 0.05, -0.02},
        frehg2::SoilParameters{0.0, 0.0, 1.0e-5, 2.0, 1.5, 0.50, 0.10, -0.02},
    };
    p.soil_id = frehg2::ReSolver::expandSurfaceSoilMap(grid, {0, 1});

    const frehg2::ReSolver solver(grid, p);
    auto state = solver.initializeLegacyState(0.0, -0.2);
    const auto clay_top = grid.getIndex(0, 0, 0);
    const auto sand_top = grid.getIndex(1, 0, 0);
    if (state.soil_id[clay_top] != 0 || state.soil_id[sand_top] != 1) {
        std::cerr << "soil ids were not assigned by column\n";
        return 1;
    }
    if (near(state.h[clay_top], state.h[sand_top], 1.0e-8)) {
        std::cerr << "different VG tables should produce different initial heads\n";
        return 1;
    }
    state.h[clay_top] = 0.1;
    state.h[sand_top] = 0.1;
    state.h[grid.getIndex(0, 0, 1)] = 0.1;
    state.h[grid.getIndex(1, 0, 1)] = 0.1;
    solver.computeConductivityFaces(state);
    if (!(state.kz[sand_top] > state.kz[clay_top])) {
        std::cerr << "higher-K soil should produce a larger vertical conductivity\n";
        return 1;
    }

    return 0;
}
