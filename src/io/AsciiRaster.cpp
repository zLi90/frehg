#include "io/AsciiRaster.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>

namespace frehg2 {

namespace {

void requireHeaderKey(const std::string& actual, const std::string& expected)
{
    if (actual != expected) {
        throw std::runtime_error("invalid ASCII raster header, expected " + expected);
    }
}

}  // namespace

void AsciiRaster::read(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open ASCII raster: " + filename);
    }

    std::string key;
    input >> key >> ncols_;
    requireHeaderKey(key, "ncols");
    input >> key >> nrows_;
    requireHeaderKey(key, "nrows");
    input >> key >> xllcorner_;
    if (key != "xllcorner" && key != "xllcenter") {
        throw std::runtime_error("invalid ASCII raster x origin header");
    }
    input >> key >> yllcorner_;
    if (key != "yllcorner" && key != "yllcenter") {
        throw std::runtime_error("invalid ASCII raster y origin header");
    }
    input >> key >> cellsize_;
    requireHeaderKey(key, "cellsize");
    input >> key >> nodata_value_;
    requireHeaderKey(key, "NODATA_value");

    if (ncols_ <= 0 || nrows_ <= 0 || cellsize_ <= 0.0) {
        throw std::runtime_error("invalid ASCII raster dimensions");
    }

    data_.assign(static_cast<std::size_t>(ncols_ * nrows_), 0.0);
    active_mask_.assign(data_.size(), 1);

    for (int row = nrows_ - 1; row >= 0; --row) {
        for (int col = 0; col < ncols_; ++col) {
            real value = 0.0;
            if (!(input >> value)) {
                throw std::runtime_error("ASCII raster ended before all cells were read");
            }
            const auto index = static_cast<std::size_t>(col + row * ncols_);
            data_[index] = value;
            active_mask_[index] = std::abs(value - nodata_value_) < 1.0e-12 ? 0 : 1;
        }
    }

#ifdef USE_KOKKOS
    data_view_ = realArr("ascii_raster_data", data_.size());
    active_mask_view_ = intArr("ascii_raster_active_mask", active_mask_.size());
    for (std::size_t i = 0; i < data_.size(); ++i) {
        data_view_(i) = data_[i];
        active_mask_view_(i) = active_mask_[i];
    }
#endif
}

int AsciiRaster::ncols() const noexcept
{
    return ncols_;
}

int AsciiRaster::nrows() const noexcept
{
    return nrows_;
}

real AsciiRaster::xllcorner() const noexcept
{
    return xllcorner_;
}

real AsciiRaster::yllcorner() const noexcept
{
    return yllcorner_;
}

real AsciiRaster::cellsize() const noexcept
{
    return cellsize_;
}

real AsciiRaster::nodataValue() const noexcept
{
    return nodata_value_;
}

index_t AsciiRaster::size() const noexcept
{
    return data_.size();
}

real AsciiRaster::value(index_t index) const
{
    if (index >= data_.size()) {
        throw std::out_of_range("ASCII raster index out of range");
    }
    return data_[index];
}

int AsciiRaster::active(index_t index) const
{
    if (index >= active_mask_.size()) {
        throw std::out_of_range("ASCII raster active-mask index out of range");
    }
    return active_mask_[index];
}

#ifdef USE_KOKKOS
const realArr& AsciiRaster::dataView() const noexcept
{
    return data_view_;
}

const intArr& AsciiRaster::activeMaskView() const noexcept
{
    return active_mask_view_;
}
#endif

}  // namespace frehg2
