// ESRI ASCII raster reader (P3.3) for DEM / bathymetry / active-mask inputs.
#ifndef FREHG2_IO_ASCII_RASTER_HPP
#define FREHG2_IO_ASCII_RASTER_HPP

#include <string>

#include "frehg2/core/define.hpp"

namespace frehg2 {

class AsciiRaster {
 public:
  // Parse an ESRI ASCII grid. Throws std::runtime_error on malformed header/data.
  static AsciiRaster fromFile(const std::string& path);
  static AsciiRaster fromString(const std::string& text);

  int ncols() const { return ncols_; }
  int nrows() const { return nrows_; }
  real xll() const { return xll_; }        // always lower-left CORNER (converted if center)
  real yll() const { return yll_; }
  real cellsize() const { return cellsize_; }
  real noData() const { return nodata_; }

  // Row-major data (j outer, i inner): index = i + j*ncols, j=0 is the first data row.
  const RealArr1DHost& data() const { return data_; }

  bool isNoData(int idx, real value) const;
  // 1.0 where valid, 0.0 where NODATA.
  RealArr1DHost activeMask() const;

 private:
  AsciiRaster() = default;

  int ncols_ = 0;
  int nrows_ = 0;
  real xll_ = 0.0;
  real yll_ = 0.0;
  real cellsize_ = 1.0;
  real nodata_ = -9999.0;
  RealArr1DHost data_;
};

}  // namespace frehg2

#endif  // FREHG2_IO_ASCII_RASTER_HPP
