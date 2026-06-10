#include "bc/Polygon.hpp"

#include <fstream>
#include <iostream>

int main()
{
    const frehg2::Polygon square({
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 2.0},
        {0.0, 2.0},
        {0.0, 0.0},
    });

    if (!square.isInside(1.0, 1.0) || !square.isInside(0.0, 1.0) ||
        square.isInside(3.0, 1.0)) {
        std::cerr << "point-in-polygon test failed\n";
        return 1;
    }

    const frehg2::Polygon disjoint({
        {3.0, 0.0},
        {4.0, 0.0},
        {4.0, 1.0},
        {3.0, 1.0},
    });
    if (square.overlaps(disjoint)) {
        std::cerr << "disjoint polygons should not overlap\n";
        return 1;
    }

    frehg2::GridSpec spec;
    spec.nx = 4;
    spec.ny = 4;
    spec.nz = 2;
    spec.dx = 1.0;
    spec.dy = 1.0;
    const frehg2::Grid grid(spec);
    if (square.selectSurfaceCells(grid).size() != 4 || square.selectGroundwaterCells(grid).size() != 8) {
        std::cerr << "polygon cell selection failed\n";
        return 1;
    }

    const std::string file = "phase9_polygon_test.poly";
    {
        std::ofstream out(file);
        out << "4\n";
        out << "0 0\n";
        out << "1 0\n";
        out << "1 1\n";
        out << "0 1\n";
    }
    const auto from_file = frehg2::Polygon::readFromFile(file);
    if (!from_file.isInside(0.5, 0.5)) {
        std::cerr << "polygon file reader failed\n";
        return 1;
    }

    return 0;
}
