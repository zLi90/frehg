#include "core/Grid.hpp"

#include <stdexcept>

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 10;
    spec.ny = 10;
    spec.nz = 3;
    spec.dx = 2.0;
    spec.dy = 3.0;
    spec.dz = 0.5;
    spec.dz_multiplier = 2.0;

    const frehg2::Grid grid(spec);

    if (grid.nSurfaceCell() != 100) {
        return 1;
    }
    if (grid.nSurfaceCellMem() != 144) {
        return 1;
    }
    if (grid.nCell() != 300) {
        return 1;
    }
    if (grid.nCellMem() != 720) {
        return 1;
    }

    const auto index = grid.getIndex(4, 5, 2);
    const auto logical = grid.getLogicalIndex(index);
    if (logical.i != 4 || logical.j != 5 || logical.k != 2) {
        return 1;
    }

    if (grid.getSurfaceIndex(0, 0) != 0) {
        return 1;
    }
    if (grid.getSurfaceIndex(1, 0) != 1) {
        return 1;
    }
    if (grid.getSurfaceIndex(0, 1) != 10) {
        return 1;
    }
    const auto surface_logical = grid.getSurfaceLogicalIndex(54);
    if (surface_logical[0] != 4 || surface_logical[1] != 5) {
        return 1;
    }

    if (grid.getIndex(0, 0, 0) != 0) {
        return 1;
    }
    if (grid.getIndex(0, 0, 1) != 1) {
        return 1;
    }
    if (grid.getIndex(1, 0, 0) != 3) {
        return 1;
    }
    if (grid.surfaceGhostIndex(frehg2::Direction::XP, 9, 3) != 100 + 20 + 3) {
        return 1;
    }
    if (grid.surfaceGhostIndex(frehg2::Direction::XM, 0, 3) != 100 + 20 + 10 + 3) {
        return 1;
    }
    if (grid.surfaceGhostIndex(frehg2::Direction::YP, 4, 9) != 100 + 4) {
        return 1;
    }
    if (grid.surfaceGhostIndex(frehg2::Direction::YM, 4, 0) != 100 + 10 + 4) {
        return 1;
    }

    if (grid.groundwaterGhostIndex(frehg2::Direction::XP, 9, 3, 2) !=
        300 + 2 * 10 * 3 + 3 * 3 + 2) {
        return 1;
    }
    if (grid.groundwaterGhostIndex(frehg2::Direction::ZP, 4, 5, 2) !=
        144 * 3 + (5 * 10 + 4)) {
        return 1;
    }
    if (grid.groundwaterGhostIndex(frehg2::Direction::ZM, 4, 5, 0) !=
        144 * 4 + (5 * 10 + 4)) {
        return 1;
    }

    if (grid.dz(0) != 0.5 || grid.dz(1) != 1.0 || grid.dz(2) != 2.0) {
        return 1;
    }

    const frehg2::Grid decomposed(spec, 2, 2, 3);
    const auto global = decomposed.localToGlobal(2, 3, 1);
    if (global.i != 12 || global.j != 13 || global.k != 1) {
        return 1;
    }
    const auto local = decomposed.globalToLocal(12, 13, 1);
    if (local.i != 2 || local.j != 3 || local.k != 1) {
        return 1;
    }

    try {
        (void)grid.getIndex(10, 0, 0);
        return 1;
    } catch (const std::out_of_range&) {
    }

    return 0;
}
