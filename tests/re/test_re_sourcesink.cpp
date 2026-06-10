#include "re/ReSolver.hpp"

#include <iostream>

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 1;
    spec.ny = 1;
    spec.nz = 8;
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
    p.bc_type = {0, 0, 0, 0, 2, 2};
    p.qbot = 1.0e-8;
    p.qtop = -1.0e-8;
    p.dt_adjust = false;

    frehg2::ReSolver solver(frehg2::Grid(spec), p);
    auto state = solver.initializeLegacyState(0.0, -0.08);
    solver.computeConductivityFaces(state);
    auto system = solver.assemblePredictorSystem(state);
    if (!(system.rhs.front() != system.rhs.back())) {
        std::cerr << "top/bottom flux source terms were not applied\n";
        return 1;
    }

    spec.nx = 3;
    spec.nz = 4;
    p.qbot = 0.0;
    p.qtop = 0.0;
    p.bc_type = {0, 0, 0, 0, 0, 2};
    p.qtop_surface = {0.0, -1.0e-8, 0.0};

    const frehg2::Grid grid(spec);
    frehg2::ReSolver polygon_solver(grid, p);
    auto polygon_state = polygon_solver.initializeLegacyState(0.0, -0.04);
    polygon_solver.computeConductivityFaces(polygon_state);
    polygon_solver.assemblePredictorSystem(polygon_state);

    const auto left_top = grid.getIndex(0, 0, 0);
    const auto middle_top = grid.getIndex(1, 0, 0);
    const auto right_top = grid.getIndex(2, 0, 0);
    if (polygon_state.grhs[left_top] != polygon_state.grhs[right_top]) {
        std::cerr << "unselected columns should retain the same zero top flux\n";
        return 1;
    }
    if (polygon_state.grhs[middle_top] == polygon_state.grhs[left_top]) {
        std::cerr << "selected polygon column did not receive distinct top flux\n";
        return 1;
    }

    frehg2::ReParameters noflow = p;
    noflow.qtop_surface = {0.0, 0.0, 0.0};
    frehg2::ReSolver noflow_solver(grid, noflow);
    auto noflow_state = noflow_solver.initializeLegacyState(0.0, -0.04);
    noflow_solver.computeConductivityFaces(noflow_state);
    noflow_solver.assemblePredictorSystem(noflow_state);
    if (noflow_state.grhs[left_top] != noflow_state.grhs[middle_top] ||
        noflow_state.gzm[left_top] != 0.0 ||
        noflow_state.gzm[middle_top] != 0.0) {
        std::cerr << "zero top fixed-flux cells should be no-flow\n";
        return 1;
    }
    return 0;
}
