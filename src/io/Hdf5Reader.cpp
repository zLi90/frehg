#include "io/Hdf5Reader.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#include <stdexcept>
#include <utility>

namespace frehg2 {

namespace {

#ifdef USE_HDF5
void checkStatus(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error(message);
    }
}
#endif

[[noreturn, maybe_unused]] void throwHdf5Unavailable()
{
    throw std::runtime_error("HDF5 support is disabled; install HDF5 under the local prefix and "
                             "configure with USE_HDF5=ON");
}

}  // namespace

Hdf5Reader::Hdf5Reader(std::string filename)
    : filename_(std::move(filename))
{
}

const std::string& Hdf5Reader::filename() const noexcept
{
    return filename_;
}

bool Hdf5Reader::isEnabled() const noexcept
{
#ifdef USE_HDF5
    return true;
#else
    return false;
#endif
}

std::vector<real> Hdf5Reader::readVector(const std::string& dataset_path) const
{
#ifdef USE_HDF5
    if (dataset_path.empty() || dataset_path.front() != '/') {
        throw std::invalid_argument("HDF5 dataset path must be absolute");
    }

    const hid_t file = H5Fopen(filename_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to open HDF5 file: " + filename_);
    }

    const hid_t dataset = H5Dopen2(file, dataset_path.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        H5Fclose(file);
        throw std::runtime_error("failed to open HDF5 dataset: " + dataset_path);
    }

    const hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0) {
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("failed to read HDF5 dataspace");
    }

    if (H5Sget_simple_extent_ndims(dataspace) != 1) {
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("HDF5 dataset is not one-dimensional: " + dataset_path);
    }

    hsize_t dims[1] = {0};
    if (H5Sget_simple_extent_dims(dataspace, dims, nullptr) < 0) {
        H5Sclose(dataspace);
        H5Dclose(dataset);
        H5Fclose(file);
        throw std::runtime_error("failed to read HDF5 dataset dimensions");
    }

    std::vector<real> values(static_cast<std::size_t>(dims[0]));
    if (!values.empty()) {
        checkStatus(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                            values.data()),
                    "failed to read HDF5 dataset");
    }
    checkStatus(H5Sclose(dataspace), "failed to close HDF5 dataspace");
    checkStatus(H5Dclose(dataset), "failed to close HDF5 dataset");
    checkStatus(H5Fclose(file), "failed to close HDF5 file");
    return values;
#else
    (void)dataset_path;
    throwHdf5Unavailable();
#endif
}

}  // namespace frehg2
