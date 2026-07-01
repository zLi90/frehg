// Polygon boundary-condition application (P12.3.3).
//
// Applies the per-region BC rule (bc_discharge / bc_depth / bc_critical) to the surface-water
// state on every owned cell tagged by the boundary PolygonIndex, at EVERY time step (the rule
// runs from Orchestrator::step via the surface-advance helper — never only at initialization).
// BCs are applied as explicit post-solve state updates, consistent with how the realized SWE
// solver applies its explicit rainfall/evaporation source (legacy evaprain runs after the KSP
// solve); the surface geometry is refreshed afterwards so the next step sees a consistent state.
//
// This unit is plain serial/host and free of PETSc types (it drives SweSolver through its public
// field accessors); on-node Kokkos/GPU is the global P9/P10 pass.
#ifndef FREHG2_BC_POLYGON_BC_HPP
#define FREHG2_BC_POLYGON_BC_HPP

#include <vector>

#include "bc/PolygonIndex.hpp"
#include "bc/PolygonRegion.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

class SweSolver;
class MpiComm;

class PolygonBC {
 public:
  // `regions` and `index` MUST be in the same polygon order (index.lookup returns the region
  // subscript). `mc` is used to compute the global per-region cell count for even distribution.
  PolygonBC(std::vector<BcRegion> regions, PolygonIndex index, const MpiComm* mc);

  bool empty() const { return regions_.empty(); }
  int nRegions() const { return static_cast<int>(regions_.size()); }

  // Apply all BC regions to the surface state for surface step `dt`. Returns the net OUTFLOW
  // volume [m^3] this step (positive = water removed from the domain; negative = injected),
  // summed over discharge/critical/depth regions. Refreshes depth + subgrid geometry when any
  // region is present. No-op (returns 0) when there are no boundary regions.
  real applySurface(SweSolver& swe, real dt) const;

 private:
  std::vector<BcRegion> regions_;
  PolygonIndex index_;
  std::vector<int> global_count_;
};

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_BC_HPP
