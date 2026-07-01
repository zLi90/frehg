// Load a 2D or 3D scalar field raster for IC application (P14.3.3).
#ifndef FREHG2_IC_RASTER_FIELD_HPP
#define FREHG2_IC_RASTER_FIELD_HPP

#include <string>
#include <vector>

#include "frehg2/core/define.hpp"

namespace frehg2 {

class Config;

// Global field in row-major order: index = gi + gj*gnx (+ gk*gnx*gny when 3D).
struct RasterField {
  int gnx = 0;
  int gny = 0;
  int gnz = 1;  // 1 => 2D column field replicated across k
  std::vector<double> values;

  double at(int gi, int gj, int gk = 0) const {
    return values[static_cast<size_t>(gi + gj * gnx + gk * gnx * gny)];
  }
};

// Resolve relative paths through Config, then load list / ESRI ASCII / HDF5 dataset.
RasterField loadRasterField(const Config& config, const std::string& file,
                            const std::string& format, const std::string& dataset, int gnx,
                            int gny, int gnz);

}  // namespace frehg2

#endif  // FREHG2_IC_RASTER_FIELD_HPP
