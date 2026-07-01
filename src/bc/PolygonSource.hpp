// Polygon source/sink application (P12.3.4).
//
// Surface sources (inflow_rate, rainfall_rate) and a subsurface sink (extraction_well), applied
// to every owned cell tagged by the source PolygonIndex at EVERY time step. Sources are explicit
// state updates applied each step (consistent with the realized SWE solver's explicit rainfall):
//   - inflow_rate    : volumetric inflow Q [m^3/s] spread evenly over the region's surface cells
//   - rainfall_rate  : extra rainfall r [m/s] over the region's surface cells (added to eta)
//   - extraction_well: Q [m^3/s] removed from the deepest GW cell of each column in the region,
//                      as an explicit water-content sink clamped at theta_r
//
// Plain serial/host; no PETSc types (drives SweSolver / ReSolver through public accessors).
#ifndef FREHG2_BC_POLYGON_SOURCE_HPP
#define FREHG2_BC_POLYGON_SOURCE_HPP

#include <vector>

#include "bc/PolygonIndex.hpp"
#include "bc/PolygonRegion.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

class SweSolver;
class ReSolver;
class MpiComm;

class PolygonSource {
 public:
  PolygonSource(std::vector<SourceRegion> regions, PolygonIndex index, const MpiComm* mc);

  bool empty() const { return regions_.empty(); }
  bool hasSurface() const { return has_surface_; }
  bool hasSubsurface() const { return has_subsurface_; }

  // Apply surface sources (inflow_rate, rainfall_rate) for surface step `dt`. Returns the net
  // water volume [m^3] added to the surface this step. Refreshes depth/geometry if applied.
  real applySurface(SweSolver& swe, real dt) const;

  // Apply subsurface extraction wells for elapsed groundwater time `dt`. Returns the net water
  // volume [m^3] removed from the subsurface this step. No-op when there are no wells.
  real applySubsurface(ReSolver& re, real dt) const;

 private:
  std::vector<SourceRegion> regions_;
  PolygonIndex index_;
  std::vector<int> global_count_;
  bool has_surface_ = false;
  bool has_subsurface_ = false;
};

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_SOURCE_HPP
