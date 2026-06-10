#include "bc/BoundaryCondition.hpp"

#include <iostream>
#include <vector>

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 3;
    spec.ny = 3;
    spec.nz = 2;
    spec.dx = 1.0;
    spec.dy = 1.0;
    const frehg2::Grid grid(spec);

    const frehg2::Polygon center({
        {1.0, 1.0},
        {2.0, 1.0},
        {2.0, 2.0},
        {1.0, 2.0},
    });

    std::vector<frehg2::real> eta(static_cast<std::size_t>(grid.nSurfaceCell()), 0.0);
    std::vector<frehg2::real> bed(static_cast<std::size_t>(grid.nSurfaceCell()), 0.25);

    const std::vector<frehg2::PolygonBoundaryCondition> surface_conditions{
        {center, frehg2::BoundaryConditionType::FIXED_WATER_LEVEL, 1.5},
    };
    frehg2::BoundaryCondition::applySurface(grid, surface_conditions, eta, bed);

    const auto center_cell = grid.getSurfaceIndex(1, 1);
    if (eta[center_cell] != 1.5) {
        std::cerr << "fixed water-level polygon BC was not applied\n";
        return 1;
    }
    if (eta[grid.getSurfaceIndex(0, 0)] != 0.0) {
        std::cerr << "polygon BC leaked outside target polygon\n";
        return 1;
    }

    const auto markers = frehg2::BoundaryCondition::markCells(grid, {
        {center, frehg2::BoundaryConditionType::FREE_OUTFLOW, 0.0},
    });
    if (markers[center_cell] == 0 || markers[grid.getSurfaceIndex(0, 0)] != 0) {
        std::cerr << "polygon BC marker assignment failed\n";
        return 1;
    }

    std::vector<frehg2::real> head(static_cast<std::size_t>(grid.nCell()), -1.0);
    frehg2::BoundaryCondition::applyGroundwaterHead(grid, {
        {center, frehg2::BoundaryConditionType::GROUNDWATER_HEAD, 0.75},
    }, head);
    if (head[grid.getIndex(1, 1, 0)] != 0.75 || head[grid.getIndex(1, 1, 1)] != 0.75) {
        std::cerr << "groundwater-head polygon BC was not applied through the column\n";
        return 1;
    }
    if (head[grid.getIndex(0, 0, 0)] != -1.0) {
        std::cerr << "groundwater-head polygon BC leaked outside target polygon\n";
        return 1;
    }

    std::vector<frehg2::real> top_flux(
        static_cast<std::size_t>(grid.nSurfaceCell()),
        0.0);
    const bool applied_flux = frehg2::BoundaryCondition::applyGroundwaterTopFlux(
        grid,
        {
            {center, frehg2::BoundaryConditionType::FIXED_FLOW_RATE, -5.0e-6},
        },
        top_flux);
    if (!applied_flux || top_flux[center_cell] != -5.0e-6) {
        std::cerr << "groundwater top-flux polygon BC was not applied\n";
        return 1;
    }
    if (top_flux[grid.getSurfaceIndex(0, 0)] != 0.0) {
        std::cerr << "groundwater top-flux polygon BC leaked outside target polygon\n";
        return 1;
    }

    return 0;
}
