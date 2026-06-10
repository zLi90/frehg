#ifndef FREHG2_IO_HDF5_READER_HPP
#define FREHG2_IO_HDF5_READER_HPP

#include "core/define.hpp"

#include <string>
#include <vector>

namespace frehg2 {

class Hdf5Reader {
public:
    explicit Hdf5Reader(std::string filename);

    const std::string& filename() const noexcept;
    bool isEnabled() const noexcept;

    std::vector<real> readVector(const std::string& dataset_path) const;

private:
    std::string filename_;
};

}  // namespace frehg2

#endif  // FREHG2_IO_HDF5_READER_HPP
