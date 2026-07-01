// Time-series reader (P3.4) for rainfall / tidal BC forcing. CSV or whitespace pairs.
#ifndef FREHG2_IO_TIME_SERIES_HPP
#define FREHG2_IO_TIME_SERIES_HPP

#include <string>
#include <vector>

#include "frehg2/core/define.hpp"

namespace frehg2 {

class TimeSeries {
 public:
  // Read "time value" pairs, one per line; columns separated by comma or whitespace.
  // Blank lines and lines beginning with '#' are ignored. Entries are sorted by time.
  static TimeSeries fromFile(const std::string& path);
  static TimeSeries fromString(const std::string& text);

  size_t size() const { return times_.size(); }
  const std::vector<real>& times() const { return times_; }
  const std::vector<real>& values() const { return values_; }

  // Linear interpolation inside the range; nearest-boundary (clamped) outside the range.
  real getValueAt(real time) const;

 private:
  TimeSeries() = default;
  std::vector<real> times_;
  std::vector<real> values_;
};

}  // namespace frehg2

#endif  // FREHG2_IO_TIME_SERIES_HPP
