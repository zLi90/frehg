#include "re/ReSolver.hpp"

#include <cmath>
#include <iostream>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

frehg2::ReParameters benchmarkParameters()
{
    frehg2::ReParameters p;
    p.ksz = 2.89e-6;
    p.ss = 1.0e-5;
    p.soil_a = 1.43;
    p.soil_n = 1.56;
    p.wcs = 0.33;
    p.wcr = 0.0;
    p.init_wc = 0.033;
    p.htop = 0.0;
    p.bc_type = {0, 0, 0, 0, 0, 1};
    return p;
}

}  // namespace

int main()
{
    const auto params = benchmarkParameters();
    const auto h = frehg2::ReSolver::headFromWaterContent(params.init_wc, params);
    if (!near(h, -42.650281, 1.0e-6)) {
        std::cerr << "unexpected VG head: " << h << "\n";
        return 1;
    }
    const auto wc = frehg2::ReSolver::waterContentFromHead(h, params);
    if (!near(wc, params.init_wc, 1.0e-12)) {
        std::cerr << "VG round-trip failed: " << wc << "\n";
        return 1;
    }

    frehg2::GridSpec spec;
    spec.nx = 1;
    spec.ny = 1;
    spec.nz = 100;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.01;
    frehg2::ReSolver solver(frehg2::Grid(spec), params);
    const auto state = solver.initializeLegacyState(0.0, -1.0);
    if (state.h.size() != 100 || !near(state.zcntr.front(), -0.005, 1.0e-12) ||
        !near(state.zcntr.back(), -0.995, 1.0e-12)) {
        std::cerr << "legacy column initialization failed\n";
        return 1;
    }
    if (!near(state.h.front(), -42.650281, 1.0e-6)) {
        std::cerr << "unexpected initial column head\n";
        return 1;
    }
    return 0;
}
