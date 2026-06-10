#include "bc/BoundaryCondition.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

}  // namespace

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 4;
    spec.ny = 2;
    spec.nz = 1;
    spec.dx = 1.0;
    spec.dy = 1.0;
    const frehg2::Grid grid(spec);

    const frehg2::Polygon left_half({
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 2.0},
        {0.0, 2.0},
    });
    const frehg2::Polygon first_cell({
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0},
    });

    std::vector<frehg2::real> depth(static_cast<std::size_t>(grid.nSurfaceCell()), 0.1);
    frehg2::BoundaryCondition::applySourceSink(grid, {
        {left_half, frehg2::SourceSinkType::RAINFALL, 0.01},
        {first_cell, frehg2::SourceSinkType::PUMPING, 0.005},
    }, depth, 10.0);

    if (!near(depth[grid.getSurfaceIndex(0, 0)], 0.15, 1.0e-12)) {
        std::cerr << "overlapping rainfall and pumping polygons produced wrong depth\n";
        return 1;
    }
    if (!near(depth[grid.getSurfaceIndex(1, 0)], 0.20, 1.0e-12) ||
        !near(depth[grid.getSurfaceIndex(1, 1)], 0.20, 1.0e-12)) {
        std::cerr << "rainfall polygon did not apply to all selected cells\n";
        return 1;
    }
    if (!near(depth[grid.getSurfaceIndex(3, 1)], 0.1, 1.0e-12)) {
        std::cerr << "source/sink leaked outside polygon\n";
        return 1;
    }

    return 0;
}
