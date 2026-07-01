// Parse polygon BC / source regions from the YAML config (P12.3.5).
//
// Frozen schema (additive; absent => no polygons => exact b1-sw / b2-gw regression):
//
//   boundaries:
//     - name: outlet
//       type: bc_discharge        # bc_discharge | bc_depth | bc_critical
//       vertices: [[x1, y1], [x2, y2], ...]   # >= 3, global coordinates
//       rate: 1.0                 # Q [m^3/s] (discharge) or h [m] (depth); unused for critical
//   sources:
//     - name: catchment_inflow
//       type: inflow_rate         # inflow_rate | rainfall_rate | extraction_well
//       vertices: [[x1, y1], ...]
//       rate: 1.0                 # Q [m^3/s] (inflow/well) or r [m/s] (rainfall)
//
// `domain.x0` / `domain.y0` (default 0) set the global origin used by PolygonIndex.
#ifndef FREHG2_BC_POLYGON_CONFIG_HPP
#define FREHG2_BC_POLYGON_CONFIG_HPP

#include <vector>

#include "bc/PolygonRegion.hpp"

namespace frehg2 {

class Config;

struct PolygonRegions {
  std::vector<BcRegion> boundaries;
  std::vector<SourceRegion> sources;
};

// Parse the `boundaries:` and `sources:` sequences. Throws std::runtime_error on a malformed
// entry (missing/unknown type, fewer than 3 vertices, non-numeric coordinate). Returns empty
// lists when the sections are absent.
PolygonRegions parsePolygonRegions(const Config& cfg);

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_CONFIG_HPP
