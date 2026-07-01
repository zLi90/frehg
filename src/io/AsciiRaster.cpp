#include "io/AsciiRaster.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace frehg2 {

namespace {
std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}
}  // namespace

AsciiRaster AsciiRaster::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("AsciiRaster: cannot open file '" + path + "'");
  }
  std::stringstream ss;
  ss << in.rdbuf();
  return fromString(ss.str());
}

AsciiRaster AsciiRaster::fromString(const std::string& text) {
  AsciiRaster r;
  std::istringstream in(text);

  bool have_ncols = false, have_nrows = false, have_x = false, have_y = false,
       have_cell = false;
  bool x_is_center = false, y_is_center = false;

  // Parse the six header lines (order-independent, case-insensitive keys).
  std::string key;
  while (!(have_ncols && have_nrows && have_x && have_y && have_cell)) {
    if (!(in >> key)) {
      throw std::runtime_error("AsciiRaster: incomplete header");
    }
    const std::string k = lower(key);
    if (k == "ncols") { in >> r.ncols_; have_ncols = true; }
    else if (k == "nrows") { in >> r.nrows_; have_nrows = true; }
    else if (k == "xllcorner") { in >> r.xll_; x_is_center = false; have_x = true; }
    else if (k == "xllcenter") { in >> r.xll_; x_is_center = true; have_x = true; }
    else if (k == "yllcorner") { in >> r.yll_; y_is_center = false; have_y = true; }
    else if (k == "yllcenter") { in >> r.yll_; y_is_center = true; have_y = true; }
    else if (k == "cellsize") { in >> r.cellsize_; have_cell = true; }
    else if (k == "nodata_value") { in >> r.nodata_; }
    else {
      throw std::runtime_error("AsciiRaster: unknown header key '" + key + "'");
    }
  }
  // An optional NODATA_value line may follow the mandatory header keys.
  // Peek: if the next token is a NODATA_value key, consume it.
  std::streampos pos = in.tellg();
  if (in >> key) {
    if (lower(key) == "nodata_value") {
      in >> r.nodata_;
    } else {
      in.clear();
      in.seekg(pos);
    }
  } else {
    in.clear();
    in.seekg(pos);
  }

  if (r.ncols_ <= 0 || r.nrows_ <= 0) {
    throw std::runtime_error("AsciiRaster: ncols/nrows must be positive");
  }
  // Convert center -> lower-left corner so xll()/yll() are consistent.
  if (x_is_center) r.xll_ -= 0.5 * r.cellsize_;
  if (y_is_center) r.yll_ -= 0.5 * r.cellsize_;

  const size_t n = static_cast<size_t>(r.ncols_) * r.nrows_;
  r.data_ = RealArr1DHost("AsciiRaster::data", n);
  for (size_t idx = 0; idx < n; ++idx) {
    double v = 0.0;
    if (!(in >> v)) {
      throw std::runtime_error("AsciiRaster: not enough data values (expected " +
                               std::to_string(n) + ")");
    }
    r.data_(idx) = static_cast<real>(v);
  }
  return r;
}

bool AsciiRaster::isNoData(int idx, real value) const {
  (void)idx;
  // Match within a small relative tolerance to be robust to text round-off.
  const real tol = std::max(real(1e-6), std::fabs(nodata_) * real(1e-9));
  return std::fabs(value - nodata_) <= tol;
}

RealArr1DHost AsciiRaster::activeMask() const {
  RealArr1DHost mask("AsciiRaster::activeMask", data_.extent(0));
  for (size_t i = 0; i < data_.extent(0); ++i) {
    mask(i) = isNoData(static_cast<int>(i), data_(i)) ? real(0.0) : real(1.0);
  }
  return mask;
}

}  // namespace frehg2
