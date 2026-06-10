#include "core/Monitor.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#include <algorithm>
#include <stdexcept>

namespace frehg2 {

namespace {

#ifdef USE_HDF5
void checkHdf5(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error(message);
    }
}
#endif

}  // namespace

real MassBalance::residual() const noexcept
{
    return new_mass - old_mass - (boundary_flux + source_sink) * dt;
}

real Monitor::pointValue(const std::vector<real>& values, index_t cell)
{
    if (cell >= values.size()) {
        throw std::out_of_range("monitor point cell is outside the value array");
    }
    return values[cell];
}

MonitorStats Monitor::summarizeCells(
    const std::vector<real>& values,
    const std::vector<index_t>& cells)
{
    if (cells.empty()) {
        return {};
    }

    MonitorStats stats;
    stats.count = cells.size();
    stats.min = pointValue(values, cells.front());
    stats.max = stats.min;
    for (const auto cell : cells) {
        const real value = pointValue(values, cell);
        stats.sum += value;
        stats.min = std::min(stats.min, value);
        stats.max = std::max(stats.max, value);
    }
    stats.mean = stats.sum / static_cast<real>(stats.count);
    return stats;
}

MassBalance Monitor::massBalance(
    real old_mass,
    real new_mass,
    real boundary_flux,
    real source_sink,
    real dt)
{
    if (dt < 0.0) {
        throw std::invalid_argument("mass-balance timestep must be non-negative");
    }
    return MassBalance{old_mass, new_mass, boundary_flux, source_sink, dt};
}

void Monitor::writePolygonHdf5(
    const std::string& filename,
    const std::string& group_name,
    const std::vector<index_t>& cells,
    const std::vector<real>& values)
{
    if (cells.size() != values.size()) {
        throw std::invalid_argument("monitor cells and values must have matching sizes");
    }

#ifdef USE_HDF5
    const hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to create monitor HDF5 file: " + filename);
    }

    const std::string group_path = group_name.empty() || group_name.front() == '/'
                                       ? group_name
                                       : "/" + group_name;
    const hid_t group = H5Gcreate2(file, group_path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) {
        H5Fclose(file);
        throw std::runtime_error("failed to create monitor HDF5 group: " + group_path);
    }

    const hsize_t dims[1] = {static_cast<hsize_t>(cells.size())};
    const hid_t space = H5Screate_simple(1, dims, nullptr);
    if (space < 0) {
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error("failed to create monitor HDF5 dataspace");
    }

    std::vector<unsigned long long> cell_values(cells.begin(), cells.end());
    hid_t dataset = H5Dcreate2(group, "cells", H5T_NATIVE_ULLONG, space, H5P_DEFAULT,
                               H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) {
        H5Sclose(space);
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error("failed to create monitor cells dataset");
    }
    if (!cell_values.empty()) {
        checkHdf5(H5Dwrite(dataset, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           cell_values.data()),
                  "failed to write monitor cells");
    }
    checkHdf5(H5Dclose(dataset), "failed to close monitor cells dataset");

    dataset = H5Dcreate2(group, "values", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT,
                         H5P_DEFAULT, H5P_DEFAULT);
    if (dataset < 0) {
        H5Sclose(space);
        H5Gclose(group);
        H5Fclose(file);
        throw std::runtime_error("failed to create monitor values dataset");
    }
    if (!values.empty()) {
        checkHdf5(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           values.data()),
                  "failed to write monitor values");
    }
    checkHdf5(H5Dclose(dataset), "failed to close monitor values dataset");
    checkHdf5(H5Sclose(space), "failed to close monitor dataspace");
    checkHdf5(H5Gclose(group), "failed to close monitor group");
    checkHdf5(H5Fclose(file), "failed to close monitor file");
#else
    (void)filename;
    (void)group_name;
    throw std::runtime_error("HDF5 support is disabled");
#endif
}

}  // namespace frehg2
