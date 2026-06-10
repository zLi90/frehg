#include "core/InitialCondition.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace frehg2 {

namespace {

struct AsciiGridData {
    int ncols = 0;
    int nrows = 0;
    int nlayers = 0;
    real nodata = -9999.0;
    std::vector<real> values;
};

std::string lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool isHeaderKey(const std::string& key)
{
    const auto normalized = lower(key);
    return normalized == "ncols" || normalized == "n_cols" || normalized == "nrows" ||
           normalized == "n_rows" || normalized == "nlayers" || normalized == "n_layers" ||
           normalized == "nz" || normalized == "xllcorner" || normalized == "xllcenter" ||
           normalized == "yllcorner" || normalized == "yllcenter" || normalized == "cellsize" ||
           normalized == "nodata_value";
}

AsciiGridData readAsciiGridData(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open ASCII IC raster: " + filename);
    }

    AsciiGridData data;
    std::string token;
    while (input >> token) {
        if (!isHeaderKey(token)) {
            data.values.push_back(static_cast<real>(std::stod(token)));
            break;
        }

        std::string raw_value;
        if (!(input >> raw_value)) {
            throw std::runtime_error("ASCII IC raster header is missing a value");
        }
        const auto key = lower(token);
        if (key == "ncols" || key == "n_cols") {
            data.ncols = std::stoi(raw_value);
        } else if (key == "nrows" || key == "n_rows") {
            data.nrows = std::stoi(raw_value);
        } else if (key == "nlayers" || key == "n_layers" || key == "nz") {
            data.nlayers = std::stoi(raw_value);
        } else if (key == "nodata_value") {
            data.nodata = static_cast<real>(std::stod(raw_value));
        }
    }

    real value = 0.0;
    while (input >> value) {
        data.values.push_back(value);
    }

    if (data.ncols <= 0) {
        throw std::runtime_error("invalid ASCII IC raster ncols header");
    }
    if (data.nrows <= 0) {
        throw std::runtime_error("invalid ASCII IC raster nrows header");
    }

    const auto cells_per_layer = static_cast<std::size_t>(data.ncols * data.nrows);
    if (cells_per_layer == 0 || data.values.empty() || data.values.size() % cells_per_layer != 0) {
        throw std::runtime_error("ASCII IC raster value count does not match raster dimensions");
    }
    const auto inferred_layers = static_cast<int>(data.values.size() / cells_per_layer);
    if (data.nlayers == 0) {
        data.nlayers = inferred_layers;
    }
    if (data.nlayers != inferred_layers) {
        throw std::runtime_error("ASCII IC raster layer count does not match value count");
    }

    return data;
}

std::size_t rawAsciiIndex(const AsciiGridData& data, int layer, int row, int col)
{
    const int file_row = data.nrows - 1 - row;
    return static_cast<std::size_t>(
        (layer * data.nrows + file_row) * data.ncols + col);
}

#ifdef USE_HDF5
void checkHdf5(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error(message);
    }
}
#endif

}  // namespace

std::vector<real> InitialCondition::loadSurfaceField(
    const Grid& grid,
    const InitialConditionSpec& spec)
{
    switch (spec.source) {
    case InitialConditionSource::CONSTANT:
        return std::vector<real>(static_cast<std::size_t>(grid.nSurfaceCell()), spec.value);
    case InitialConditionSource::ASCII_RASTER:
        return readAsciiSurface(grid, spec.filename);
    case InitialConditionSource::HDF5_VECTOR: {
        auto values = readHdf5Vector(spec.filename, spec.dataset);
        if (values.size() != static_cast<std::size_t>(grid.nSurfaceCell())) {
            throw std::runtime_error("surface HDF5 IC size does not match the grid");
        }
        return values;
    }
    }
    throw std::invalid_argument("unknown initial condition source");
}

std::vector<real> InitialCondition::loadGroundwaterField(
    const Grid& grid,
    const InitialConditionSpec& spec)
{
    if (spec.source == InitialConditionSource::CONSTANT) {
        return std::vector<real>(static_cast<std::size_t>(grid.nCell()), spec.value);
    }

    if (spec.source == InitialConditionSource::ASCII_RASTER) {
        return readAsciiGroundwaterField(grid, spec.filename);
    }

    auto values = readHdf5Vector(spec.filename, spec.dataset);
    if (values.size() != static_cast<std::size_t>(grid.nCell())) {
        throw std::runtime_error("groundwater HDF5 IC size does not match the grid");
    }
    return values;
}

std::vector<real> InitialCondition::readAsciiGroundwaterField(
    const Grid& grid,
    const std::string& filename)
{
    const auto data = readAsciiGridData(filename);
    if (data.ncols != grid.nx() || data.nrows != grid.ny()) {
        throw std::runtime_error("ASCII groundwater raster dimensions do not match surface grid");
    }

    std::vector<real> values(static_cast<std::size_t>(grid.nCell()), 0.0);
    if (data.nlayers == 1) {
        for (int j = 0; j < grid.ny(); ++j) {
            for (int i = 0; i < grid.nx(); ++i) {
                const auto value = data.values[rawAsciiIndex(data, 0, j, i)];
                if (std::abs(value - data.nodata) < 1.0e-12) {
                    throw std::runtime_error("ASCII groundwater raster contains NODATA");
                }
                for (int k = 0; k < grid.nz(); ++k) {
                    values[grid.getIndex(i, j, k)] = value;
                }
            }
        }
        return values;
    }

    if (data.nlayers != grid.nz()) {
        throw std::runtime_error("ASCII groundwater raster layer count does not match grid nz");
    }
    for (int k = 0; k < grid.nz(); ++k) {
        for (int j = 0; j < grid.ny(); ++j) {
            for (int i = 0; i < grid.nx(); ++i) {
                const auto value = data.values[rawAsciiIndex(data, k, j, i)];
                if (std::abs(value - data.nodata) < 1.0e-12) {
                    throw std::runtime_error("ASCII groundwater raster contains NODATA");
                }
                values[grid.getIndex(i, j, k)] = value;
            }
        }
    }
    return values;
}

std::vector<real> InitialCondition::readAsciiSurface(
    const Grid& grid,
    const std::string& filename)
{
    const auto data = readAsciiGridData(filename);
    if (data.ncols != grid.nx() || data.nrows != grid.ny() || data.nlayers != 1) {
        throw std::runtime_error("ASCII IC raster dimensions do not match surface grid");
    }

    std::vector<real> values(static_cast<std::size_t>(grid.nSurfaceCell()), 0.0);
    for (int row = 0; row < grid.ny(); ++row) {
        for (int col = 0; col < grid.nx(); ++col) {
            const auto value = data.values[rawAsciiIndex(data, 0, row, col)];
            if (std::abs(value - data.nodata) < 1.0e-12) {
                throw std::runtime_error("ASCII IC raster contains NODATA in an active cell");
            }
            values[grid.getSurfaceIndex(col, row)] = value;
        }
    }
    return values;
}

std::vector<real> InitialCondition::readHdf5Vector(
    const std::string& filename,
    const std::string& dataset)
{
#ifdef USE_HDF5
    const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to open HDF5 IC file: " + filename);
    }
    const hid_t dset = H5Dopen2(file, dataset.c_str(), H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file);
        throw std::runtime_error("failed to open HDF5 IC dataset: " + dataset);
    }
    const hid_t space = H5Dget_space(dset);
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 1) {
        H5Dclose(dset);
        H5Fclose(file);
        throw std::runtime_error("HDF5 IC dataset must be one-dimensional");
    }
    hsize_t dims[1] = {0};
    checkHdf5(H5Sget_simple_extent_dims(space, dims, nullptr),
              "failed to read HDF5 IC dimensions");
    std::vector<real> values(static_cast<std::size_t>(dims[0]));
    if (!values.empty()) {
        checkHdf5(H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          values.data()),
                  "failed to read HDF5 IC values");
    }
    checkHdf5(H5Sclose(space), "failed to close HDF5 IC dataspace");
    checkHdf5(H5Dclose(dset), "failed to close HDF5 IC dataset");
    checkHdf5(H5Fclose(file), "failed to close HDF5 IC file");
    return values;
#else
    (void)filename;
    (void)dataset;
    throw std::runtime_error("HDF5 support is disabled");
#endif
}

}  // namespace frehg2
