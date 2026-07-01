#include "io/TimeSeries.hpp"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace frehg2 {

TimeSeries TimeSeries::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("TimeSeries: cannot open file '" + path + "'");
  }
  std::stringstream ss;
  ss << in.rdbuf();
  return fromString(ss.str());
}

TimeSeries TimeSeries::fromString(const std::string& text) {
  TimeSeries ts;
  std::istringstream in(text);
  std::string line;
  std::vector<std::pair<real, real>> rows;
  while (std::getline(in, line)) {
    // Strip a trailing comment and skip blanks.
    const auto hash = line.find('#');
    if (hash != std::string::npos) line = line.substr(0, hash);
    for (char& c : line) {
      if (c == ',') c = ' ';
    }
    std::istringstream ls(line);
    double t = 0.0, v = 0.0;
    if (ls >> t >> v) {
      rows.emplace_back(static_cast<real>(t), static_cast<real>(v));
    }
  }
  if (rows.empty()) {
    throw std::runtime_error("TimeSeries: no valid (time value) rows parsed");
  }
  std::stable_sort(rows.begin(), rows.end(),
                   [](const auto& a, const auto& b) { return a.first < b.first; });
  ts.times_.reserve(rows.size());
  ts.values_.reserve(rows.size());
  for (const auto& r : rows) {
    ts.times_.push_back(r.first);
    ts.values_.push_back(r.second);
  }
  return ts;
}

real TimeSeries::getValueAt(real time) const {
  if (times_.empty()) {
    throw std::runtime_error("TimeSeries: empty series");
  }
  if (time <= times_.front()) return values_.front();  // clamp before start
  if (time >= times_.back()) return values_.back();     // clamp after end

  // First sample strictly greater than `time`.
  const auto it = std::upper_bound(times_.begin(), times_.end(), time);
  const size_t hi = static_cast<size_t>(it - times_.begin());
  const size_t lo = hi - 1;
  const real t0 = times_[lo], t1 = times_[hi];
  const real v0 = values_[lo], v1 = values_[hi];
  if (t1 == t0) return v0;
  const real w = (time - t0) / (t1 - t0);
  return v0 + w * (v1 - v0);
}

}  // namespace frehg2
