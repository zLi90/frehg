#include "re/ReSolver.hpp"

#include <cmath>
#include <iostream>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

frehg2::ReSolver makeSolver()
{
    frehg2::GridSpec spec;
    spec.nx = 1;
    spec.ny = 1;
    spec.nz = 10;
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
    p.bc_type = {0, 0, 0, 0, 0, 1};
    return frehg2::ReSolver(frehg2::Grid(spec), p);
}

}  // namespace

int main()
{
    auto solver = makeSolver();
    auto state = solver.initializeLegacyState(0.0, -0.1);
    state.h[0] = -1.0;
    state.h[1] = -2.0;
    solver.computeConductivityFaces(state);
    const auto expected_face = frehg2::ReSolver::conductivityFace(
        state.h[0], state.h[1], solver.parameters().ksz, solver.parameters());
    if (!near(state.kz[0], expected_face, 1.0e-18)) {
        std::cerr << "arithmetic K-face mismatch\n";
        return 1;
    }

    const auto system = solver.assemblePredictorSystem(state);
    if (system.n != 10 || system.entries.empty()) {
        std::cerr << "empty RE matrix\n";
        return 1;
    }
    bool has_top_diag = false;
    bool has_top_lower = false;
    for (const auto& entry : system.entries) {
        if (entry.row == 0 && entry.col == 0) {
            has_top_diag = true;
        }
        if (entry.row == 0 && entry.col == 1) {
            has_top_lower = true;
        }
    }
    if (!has_top_diag || !has_top_lower) {
        std::cerr << "top row matrix structure mismatch\n";
        return 1;
    }
    if (!(state.gct[0] > 0.0 && state.gzp[0] < 0.0 && state.gzm[0] < 0.0)) {
        std::cerr << "unexpected RE coefficients\n";
        return 1;
    }

    frehg2::GridSpec spec;
    spec.nx = 2;
    spec.ny = 1;
    spec.nz = 3;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.1;

    frehg2::ReParameters p = solver.parameters();
    p.ksx = p.ksz;
    p.ksy = p.ksz;
    p.bc_type = {0, 0, 0, 0, 0, 0};
    frehg2::ReSolver full3d_solver(frehg2::Grid(spec), p);
    auto full3d_state = full3d_solver.initializeLegacyState(0.0, -0.3);
    full3d_solver.computeConductivityFaces(full3d_state);
    const auto full3d_system = full3d_solver.assemblePredictorSystem(full3d_state);
    const int west_top = static_cast<int>(full3d_solver.grid().getIndex(0, 0, 0));
    const int east_top = static_cast<int>(full3d_solver.grid().getIndex(1, 0, 0));
    bool has_x_coupling = false;
    for (const auto& entry : full3d_system.entries) {
        if (entry.row == west_top && entry.col == east_top && entry.value < 0.0) {
            has_x_coupling = true;
        }
    }
    if (!has_x_coupling || !(full3d_state.kx[west_top] > 0.0)) {
        std::cerr << "full 3D RE matrix is missing x-direction coupling\n";
        return 1;
    }
    return 0;
}
