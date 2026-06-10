#include "io/Hdf5Writer.hpp"

#ifdef USE_HDF5
#include <hdf5.h>
#endif

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace frehg2 {

namespace {

#ifdef USE_HDF5
bool linkExists(hid_t location, const std::string& path)
{
    H5E_auto2_t old_func = nullptr;
    void* old_client_data = nullptr;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    const htri_t exists = H5Lexists(location, path.c_str(), H5P_DEFAULT);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
    return exists > 0;
}

void checkStatus(herr_t status, const std::string& message)
{
    if (status < 0) {
        throw std::runtime_error(message);
    }
}

void ensureGroups(hid_t file, const std::string& dataset_path)
{
    std::size_t start = 0;
    while (true) {
        const auto slash = dataset_path.find('/', start + 1);
        if (slash == std::string::npos) {
            break;
        }

        const auto group_path = dataset_path.substr(0, slash);
        if (!group_path.empty() && !linkExists(file, group_path)) {
            const hid_t group = H5Gcreate2(
                file, group_path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (group < 0) {
                throw std::runtime_error("failed to create HDF5 group: " + group_path);
            }
            checkStatus(H5Gclose(group), "failed to close HDF5 group");
        }
        start = slash;
    }
}

void writeStringAttribute(hid_t object, const std::string& name, const std::string& value)
{
    const hid_t string_type = H5Tcopy(H5T_C_S1);
    if (string_type < 0) {
        throw std::runtime_error("failed to create HDF5 string type");
    }
    checkStatus(H5Tset_size(string_type, H5T_VARIABLE), "failed to set HDF5 string size");

    const hid_t scalar_space = H5Screate(H5S_SCALAR);
    if (scalar_space < 0) {
        H5Tclose(string_type);
        throw std::runtime_error("failed to create HDF5 scalar dataspace");
    }

    if (linkExists(object, name)) {
        checkStatus(H5Adelete(object, name.c_str()), "failed to delete HDF5 attribute");
    }

    const hid_t attribute = H5Acreate2(
        object, name.c_str(), string_type, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
    if (attribute < 0) {
        H5Sclose(scalar_space);
        H5Tclose(string_type);
        throw std::runtime_error("failed to create HDF5 attribute: " + name);
    }

    const char* raw_value = value.c_str();
    checkStatus(H5Awrite(attribute, string_type, &raw_value), "failed to write HDF5 attribute");
    checkStatus(H5Aclose(attribute), "failed to close HDF5 attribute");
    checkStatus(H5Sclose(scalar_space), "failed to close HDF5 dataspace");
    checkStatus(H5Tclose(string_type), "failed to close HDF5 string type");
}
#endif

[[noreturn, maybe_unused]] void throwHdf5Unavailable()
{
    throw std::runtime_error("HDF5 support is disabled; install HDF5 under the local prefix and "
                             "configure with USE_HDF5=ON");
}

}  // namespace

Hdf5Writer::Hdf5Writer(std::string filename)
    : filename_(std::move(filename))
{
#ifdef USE_HDF5
    const hid_t file = H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to create HDF5 file: " + filename_);
    }
    checkStatus(H5Fclose(file), "failed to close HDF5 file");
#endif
}

const std::string& Hdf5Writer::filename() const noexcept
{
    return filename_;
}

bool Hdf5Writer::isEnabled() const noexcept
{
#ifdef USE_HDF5
    return true;
#else
    return false;
#endif
}

void Hdf5Writer::writeMetadata(const std::string& title, const std::string& version)
{
#ifdef USE_HDF5
    const hid_t file = H5Fopen(filename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to open HDF5 file: " + filename_);
    }

    const hid_t simulation_group =
        linkExists(file, "/simulation")
            ? H5Gopen2(file, "/simulation", H5P_DEFAULT)
            : H5Gcreate2(file, "/simulation", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (simulation_group < 0) {
        H5Fclose(file);
        throw std::runtime_error("failed to open or create /simulation group");
    }

    writeStringAttribute(simulation_group, "title", title);
    writeStringAttribute(simulation_group, "version", version);
    checkStatus(H5Gclose(simulation_group), "failed to close /simulation group");
    checkStatus(H5Fclose(file), "failed to close HDF5 file");
#else
    (void)title;
    (void)version;
    throwHdf5Unavailable();
#endif
}

void Hdf5Writer::writeVector(const std::string& dataset_path, const std::vector<real>& values)
{
#ifdef USE_HDF5
    if (dataset_path.empty() || dataset_path.front() != '/') {
        throw std::invalid_argument("HDF5 dataset path must be absolute");
    }

    const hid_t file = H5Fopen(filename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to open HDF5 file: " + filename_);
    }
    ensureGroups(file, dataset_path);

    if (linkExists(file, dataset_path)) {
        checkStatus(H5Ldelete(file, dataset_path.c_str(), H5P_DEFAULT),
                    "failed to delete existing HDF5 dataset");
    }

    const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    const hid_t dataspace = H5Screate_simple(1, dims, nullptr);
    if (dataspace < 0) {
        H5Fclose(file);
        throw std::runtime_error("failed to create HDF5 dataspace");
    }

    hid_t properties = H5P_DEFAULT;
    if (!values.empty()) {
        properties = H5Pcreate(H5P_DATASET_CREATE);
        if (properties < 0) {
            H5Sclose(dataspace);
            H5Fclose(file);
            throw std::runtime_error("failed to create HDF5 dataset properties");
        }
        const hsize_t chunk_dims[1] = {std::min<hsize_t>(dims[0], 1024)};
        checkStatus(H5Pset_chunk(properties, 1, chunk_dims), "failed to set HDF5 chunks");
        checkStatus(H5Pset_deflate(properties, 6), "failed to set HDF5 compression");
    }

    const hid_t dataset = H5Dcreate2(
        file, dataset_path.c_str(), H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, properties,
        H5P_DEFAULT);
    if (dataset < 0) {
        if (properties != H5P_DEFAULT) {
            H5Pclose(properties);
        }
        H5Sclose(dataspace);
        H5Fclose(file);
        throw std::runtime_error("failed to create HDF5 dataset: " + dataset_path);
    }

    if (!values.empty()) {
        checkStatus(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                             values.data()),
                    "failed to write HDF5 dataset");
    }
    checkStatus(H5Dclose(dataset), "failed to close HDF5 dataset");
    if (properties != H5P_DEFAULT) {
        checkStatus(H5Pclose(properties), "failed to close HDF5 properties");
    }
    checkStatus(H5Sclose(dataspace), "failed to close HDF5 dataspace");
    checkStatus(H5Fclose(file), "failed to close HDF5 file");
#else
    (void)dataset_path;
    (void)values;
    throwHdf5Unavailable();
#endif
}

void Hdf5Writer::writeIntVector(const std::string& dataset_path, const std::vector<int>& values)
{
#ifdef USE_HDF5
    if (dataset_path.empty() || dataset_path.front() != '/') {
        throw std::invalid_argument("HDF5 dataset path must be absolute");
    }

    const hid_t file = H5Fopen(filename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file < 0) {
        throw std::runtime_error("failed to open HDF5 file: " + filename_);
    }
    ensureGroups(file, dataset_path);

    if (linkExists(file, dataset_path)) {
        checkStatus(H5Ldelete(file, dataset_path.c_str(), H5P_DEFAULT),
                    "failed to delete existing HDF5 dataset");
    }

    const hsize_t dims[1] = {static_cast<hsize_t>(values.size())};
    const hid_t dataspace = H5Screate_simple(1, dims, nullptr);
    if (dataspace < 0) {
        H5Fclose(file);
        throw std::runtime_error("failed to create HDF5 dataspace");
    }

    const hid_t dataset = H5Dcreate2(
        file, dataset_path.c_str(), H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT,
        H5P_DEFAULT);
    if (dataset < 0) {
        H5Sclose(dataspace);
        H5Fclose(file);
        throw std::runtime_error("failed to create HDF5 dataset: " + dataset_path);
    }
    if (!values.empty()) {
        checkStatus(H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                             values.data()),
                    "failed to write HDF5 dataset");
    }
    checkStatus(H5Dclose(dataset), "failed to close HDF5 dataset");
    checkStatus(H5Sclose(dataspace), "failed to close HDF5 dataspace");
    checkStatus(H5Fclose(file), "failed to close HDF5 file");
#else
    (void)dataset_path;
    (void)values;
    throwHdf5Unavailable();
#endif
}

}  // namespace frehg2
