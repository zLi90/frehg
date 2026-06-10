#include "bc/Polygon.hpp"
#include "core/Monitor.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

#ifdef USE_HDF5
std::vector<frehg2::real> readValues(const std::string& filename)
{
    const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    const hid_t dataset = H5Dopen2(file, "/monitor/values", H5P_DEFAULT);
    const hid_t space = H5Dget_space(dataset);
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(space, dims, nullptr);
    std::vector<frehg2::real> values(static_cast<std::size_t>(dims[0]));
    H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
    return values;
}
#endif

}  // namespace

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 3;
    spec.ny = 2;
    spec.nz = 1;
    spec.dx = 1.0;
    spec.dy = 1.0;
    const frehg2::Grid grid(spec);

    const std::vector<frehg2::real> values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    if (!near(frehg2::Monitor::pointValue(values, grid.getSurfaceIndex(1, 0)), 2.0, 0.0)) {
        std::cerr << "point monitor failed\n";
        return 1;
    }

    const frehg2::Polygon left_half({
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 2.0},
        {0.0, 2.0},
    });
    const auto cells = left_half.selectSurfaceCells(grid);
    const auto stats = frehg2::Monitor::summarizeCells(values, cells);
    if (stats.count != 4 || !near(stats.sum, 12.0, 1.0e-12) ||
        !near(stats.mean, 3.0, 1.0e-12) || !near(stats.min, 1.0, 1.0e-12) ||
        !near(stats.max, 5.0, 1.0e-12)) {
        std::cerr << "polygon monitor statistics failed\n";
        return 1;
    }

    const auto balance = frehg2::Monitor::massBalance(10.0, 13.0, 0.2, 0.1, 10.0);
    if (!near(balance.residual(), 0.0, 1.0e-12)) {
        std::cerr << "mass-balance residual failed\n";
        return 1;
    }

#ifdef USE_HDF5
    const std::string filename = "phase11_monitor_polygon.h5";
    std::vector<frehg2::real> selected_values;
    for (const auto cell : cells) {
        selected_values.push_back(values[cell]);
    }
    frehg2::Monitor::writePolygonHdf5(filename, "monitor", cells, selected_values);
    if (readValues(filename) != selected_values) {
        std::cerr << "polygon monitor HDF5 output failed\n";
        return 1;
    }
#endif

    return 0;
}
