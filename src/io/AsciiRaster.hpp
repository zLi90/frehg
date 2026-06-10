#ifndef FREHG2_IO_ASCII_RASTER_HPP
#define FREHG2_IO_ASCII_RASTER_HPP

#include "core/define.hpp"

#include <string>
#include <vector>

namespace frehg2 {

class AsciiRaster {
public:
    void read(const std::string& filename);

    int ncols() const noexcept;
    int nrows() const noexcept;
    real xllcorner() const noexcept;
    real yllcorner() const noexcept;
    real cellsize() const noexcept;
    real nodataValue() const noexcept;

    index_t size() const noexcept;
    real value(index_t index) const;
    int active(index_t index) const;

#ifdef USE_KOKKOS
    const realArr& dataView() const noexcept;
    const intArr& activeMaskView() const noexcept;
#endif

private:
    int ncols_ = 0;
    int nrows_ = 0;
    real xllcorner_ = 0.0;
    real yllcorner_ = 0.0;
    real cellsize_ = 0.0;
    real nodata_value_ = -9999.0;

    std::vector<real> data_;
    std::vector<int> active_mask_;

#ifdef USE_KOKKOS
    realArr data_view_;
    intArr active_mask_view_;
#endif
};

}  // namespace frehg2

#endif  // FREHG2_IO_ASCII_RASTER_HPP
