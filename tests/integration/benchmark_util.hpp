// Shared helpers for the P18 b3-b6 benchmark integration tests. Each test loads a committed v2
// benchmark YAML, caps the step count (the real runs are far too long for CI), redirects output
// into the test scratch dir, and asserts review-tier physics (the registry gates b3-b6 `review`,
// not strict legacy L2).
#pragma once

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

#include "core/Grid.hpp"
#include "io/Config.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {
namespace bench {

// Load a benchmark YAML, cap time.max_steps, and point output at `out_h5`. The returned Config
// keeps `config_dir` = the benchmark directory so relative DEM/rain/bath paths still resolve.
inline Config loadCapped(const std::string& yaml_path, const std::string& config_dir,
                         long long max_steps, const std::string& out_h5) {
  YAML::Node n = YAML::LoadFile(yaml_path);
  n["time"]["max_steps"] = max_steps;
  if (!n["output"]) n["output"] = YAML::Node(YAML::NodeType::Map);
  n["output"]["filename"] = out_h5;
  std::stringstream ss;
  ss << n;
  return Config::fromString(ss.str(), config_dir);
}

// Total ponded surface-water volume [m^3] over owned cells (depth = max(eta-bottom,0)).
inline double surfaceWaterVolume(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const auto& f = swe.fields();
  const double cellA = g.dx() * g.dy();
  double vol = 0.0;
  bool finite = true;
  for (int j = 0; j < g.ny(); ++j)
    for (int i = 0; i < g.nx(); ++i) {
      const int c = g.getSurfaceIndex(i, j);
      const double depth = std::max(static_cast<double>(f.eta(c) - f.bottom(c)), 0.0);
      if (!std::isfinite(depth)) finite = false;
      vol += depth * cellA;
    }
  return finite ? vol : std::nan("");
}

// Max ponded depth over owned cells (NaN if any non-finite).
inline double maxSurfaceDepth(const SweSolver& swe) {
  const Grid& g = swe.grid();
  const auto& f = swe.fields();
  double m = 0.0;
  for (int j = 0; j < g.ny(); ++j)
    for (int i = 0; i < g.nx(); ++i) {
      const int c = g.getSurfaceIndex(i, j);
      const double d = static_cast<double>(f.eta(c) - f.bottom(c));
      if (!std::isfinite(d)) return std::nan("");
      m = std::max(m, d);
    }
  return m;
}

}  // namespace bench
}  // namespace frehg2
