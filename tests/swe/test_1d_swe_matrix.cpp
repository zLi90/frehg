#include "swe/SweSolver.hpp"

#include <cmath>
#include <vector>

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 10;
    spec.ny = 1;
    spec.nz = 1;
    spec.dx = 2.0;
    spec.dy = 1.0;

    frehg2::SweParameters params;
    params.dt = 0.5;
    params.gravity = 9.81;

    const frehg2::Grid grid(spec);
    const frehg2::SweSolver solver(grid, params);

    const auto nmem = static_cast<std::size_t>(grid.nSurfaceCellMem());
    std::vector<frehg2::real> eta(nmem, 1.0);
    std::vector<frehg2::real> bed(nmem, 0.0);
    std::vector<int> active(nmem, 1);
    std::vector<frehg2::real> ex(nmem, 0.0);
    std::vector<frehg2::real> ey(nmem, 0.0);
    std::vector<frehg2::real> dx_factor(nmem, 1.0);
    std::vector<frehg2::real> dy_factor(nmem, 0.0);

    const auto system = solver.assembleLinearSystem(eta, bed, active, ex, ey, dx_factor, dy_factor);
    if (system.n != 10) {
        return 1;
    }

    const frehg2::real area = spec.dx * spec.dy;
    const frehg2::real offdiag = params.gravity * params.dt * params.dt * area;
    const frehg2::real interior_diag = area + 2.0 * offdiag;
    const frehg2::real boundary_diag = area + offdiag;

    int diagonal_count = 0;
    int offdiag_count = 0;
    for (const auto& entry : system.entries) {
        if (entry.row == entry.col) {
            ++diagonal_count;
            const auto expected = (entry.row == 0 || entry.row == 9) ? boundary_diag : interior_diag;
            if (std::abs(entry.value - expected) > 1.0e-12) {
                return 1;
            }
        } else {
            ++offdiag_count;
            if (std::abs(entry.value + offdiag) > 1.0e-12) {
                return 1;
            }
            if (std::abs(entry.row - entry.col) != 1) {
                return 1;
            }
        }
    }

    if (diagonal_count != 10 || offdiag_count != 18) {
        return 1;
    }
    for (const auto rhs : system.rhs) {
        if (std::abs(rhs - area) > 1.0e-12) {
            return 1;
        }
    }

    return 0;
}
