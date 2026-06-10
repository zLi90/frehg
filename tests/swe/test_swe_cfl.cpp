#include "swe/SweSolver.hpp"

#include <cmath>
#include <vector>

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 2;
    spec.ny = 1;
    spec.nz = 1;
    spec.dx = 4.0;
    spec.dy = 4.0;

    frehg2::SweParameters params;
    params.dt = 0.5;
    params.gravity = 9.0;

    const frehg2::SweSolver solver(frehg2::Grid(spec), params);
    const std::vector<frehg2::real> h = {1.0, 1.0};
    const std::vector<frehg2::real> u = {1.0, 0.0};
    const std::vector<frehg2::real> v = {0.0, 0.0};

    const auto cfl = solver.computeCflDiagnostic(h, u, v);
    const auto expected = 0.5 * (1.0 + 3.0) / 4.0;
    return std::abs(cfl - expected) < 1.0e-12 ? 0 : 1;
}
