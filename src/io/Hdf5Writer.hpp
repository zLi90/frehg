#ifndef FREHG2_IO_HDF5_WRITER_HPP
#define FREHG2_IO_HDF5_WRITER_HPP

#include "core/define.hpp"

#include <string>
#include <vector>

namespace frehg2 {

class Hdf5Writer {
public:
    explicit Hdf5Writer(std::string filename);

    const std::string& filename() const noexcept;
    bool isEnabled() const noexcept;

    void writeMetadata(const std::string& title, const std::string& version);
    void writeVector(const std::string& dataset_path, const std::vector<real>& values);
    void writeIntVector(const std::string& dataset_path, const std::vector<int>& values);

private:
    std::string filename_;
};

}  // namespace frehg2

#endif  // FREHG2_IO_HDF5_WRITER_HPP
