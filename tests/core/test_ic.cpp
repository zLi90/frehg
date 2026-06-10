#include "core/InitialCondition.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

#ifdef USE_HDF5
void writeVectorHdf5(const std::string& filename, const std::vector<frehg2::real>& values)
{
    const hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    const hid_t space = H5Screate_simple(1, dims, nullptr);
    const hid_t dataset =
        H5Dcreate2(file, "/ic", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
    H5Dclose(dataset);
    H5Sclose(space);
    H5Fclose(file);
}
#endif

}  // namespace

int main()
{
    frehg2::GridSpec spec;
    spec.nx = 2;
    spec.ny = 2;
    spec.nz = 3;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.1;
    const frehg2::Grid grid(spec);

    const auto constant = frehg2::InitialCondition::loadSurfaceField(
        grid,
        frehg2::InitialConditionSpec{
            "h", frehg2::InitialConditionSource::CONSTANT, 2.5, "", ""});
    if (constant.size() != 4 || !near(constant[0], 2.5, 0.0) ||
        !near(constant[3], 2.5, 0.0)) {
        std::cerr << "constant surface IC failed\n";
        return 1;
    }

    const std::string ascii_file = "phase11_ic.asc";
    {
        std::ofstream out(ascii_file);
        out << "ncols 2\n";
        out << "nrows 2\n";
        out << "xllcorner 0\n";
        out << "yllcorner 0\n";
        out << "cellsize 1\n";
        out << "NODATA_value -9999\n";
        out << "3 4\n";
        out << "1 2\n";
    }
    const auto raster = frehg2::InitialCondition::loadSurfaceField(
        grid,
        frehg2::InitialConditionSpec{
            "h", frehg2::InitialConditionSource::ASCII_RASTER, 0.0, ascii_file, ""});
    if (!near(raster[grid.getSurfaceIndex(0, 0)], 1.0, 0.0) ||
        !near(raster[grid.getSurfaceIndex(1, 1)], 4.0, 0.0)) {
        std::cerr << "ASCII raster IC orientation failed\n";
        return 1;
    }

    const auto gw_from_raster = frehg2::InitialCondition::loadGroundwaterField(
        grid,
        frehg2::InitialConditionSpec{
            "wc", frehg2::InitialConditionSource::ASCII_RASTER, 0.0, ascii_file, ""});
    if (gw_from_raster.size() != static_cast<std::size_t>(grid.nCell()) ||
        !near(gw_from_raster[grid.getIndex(1, 1, 2)], 4.0, 0.0)) {
        std::cerr << "groundwater raster IC expansion failed\n";
        return 1;
    }

    const std::string ascii_3d_file = "phase11_ic_3d.asc";
    {
        std::ofstream out(ascii_3d_file);
        out << "ncols 2\n";
        out << "nrows 2\n";
        out << "nlayers 3\n";
        out << "xllcorner 0\n";
        out << "yllcorner 0\n";
        out << "cellsize 1\n";
        out << "NODATA_value -9999\n";
        out << "3 4\n";
        out << "1 2\n";
        out << "7 8\n";
        out << "5 6\n";
        out << "11 12\n";
        out << "9 10\n";
    }
    const auto gw_from_3d_raster = frehg2::InitialCondition::loadGroundwaterField(
        grid,
        frehg2::InitialConditionSpec{
            "h", frehg2::InitialConditionSource::ASCII_RASTER, 0.0, ascii_3d_file, ""});
    if (!near(gw_from_3d_raster[grid.getIndex(0, 0, 0)], 1.0, 0.0) ||
        !near(gw_from_3d_raster[grid.getIndex(1, 1, 1)], 8.0, 0.0) ||
        !near(gw_from_3d_raster[grid.getIndex(0, 0, 2)], 9.0, 0.0)) {
        std::cerr << "3D groundwater ASCII raster IC orientation failed\n";
        return 1;
    }

#ifdef USE_HDF5
    const std::string h5_file = "phase11_ic.h5";
    writeVectorHdf5(h5_file, constant);
    const auto from_hdf5 = frehg2::InitialCondition::loadSurfaceField(
        grid,
        frehg2::InitialConditionSpec{
            "h", frehg2::InitialConditionSource::HDF5_VECTOR, 0.0, h5_file, "/ic"});
    if (from_hdf5 != constant) {
        std::cerr << "HDF5 IC vector load failed\n";
        return 1;
    }
#endif

    return 0;
}
