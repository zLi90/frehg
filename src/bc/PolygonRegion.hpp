// Typed polygon regions for boundary conditions and source/sink terms (P12).
//
// A region is a Polygon plus a rule (its `kind`) and a scalar `rate`. Boundary regions live in
// the YAML `boundaries:` list; source/sink regions live in the `sources:` list. The polygon
// `bc_type` integer arrays of legacy `frehg` are DEPRECATED (will be removed in P23) — this is
// the primary BC path from P12 onward.
#ifndef FREHG2_BC_POLYGON_REGION_HPP
#define FREHG2_BC_POLYGON_REGION_HPP

#include <string>

#include "bc/Polygon.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

// Boundary-condition rules applied to surface cells inside the polygon.
enum class BcKind {
  Discharge,  // bc_discharge: prescribe net volumetric flux Q [m^3/s] (>0 = outflow)
  Depth,      // bc_depth: prescribe ponded depth h [m] (overrides the solved eta)
  Critical    // bc_critical: critical-depth weir outflow q = sqrt(g) h^{3/2} per unit width
};

// Source/sink rules.
enum class SourceKind {
  InflowRate,      // inflow_rate: volumetric inflow Q [m^3/s] spread evenly over surface cells
  RainfallRate,    // rainfall_rate: extra rainfall r [m/s] over surface cells in the polygon
  ExtractionWell   // extraction_well: Q [m^3/s] (>0 = extraction) from the deepest GW cell/column
};

struct BcRegion {
  Polygon polygon;
  BcKind kind = BcKind::Discharge;
  real rate = 0.0;  // Q [m^3/s] for Discharge; h [m] for Depth; unused for Critical
};

struct SourceRegion {
  Polygon polygon;
  SourceKind kind = SourceKind::InflowRate;
  real rate = 0.0;  // Q [m^3/s] (InflowRate / ExtractionWell) or r [m/s] (RainfallRate)
};

// String <-> enum helpers (throw std::runtime_error on an unknown type string).
BcKind bcKindFromString(const std::string& s);
SourceKind sourceKindFromString(const std::string& s);
const char* toString(BcKind k);
const char* toString(SourceKind k);

// Whether a YAML `type:` string names a boundary rule or a source rule.
bool isBcTypeString(const std::string& s);
bool isSourceTypeString(const std::string& s);

}  // namespace frehg2

#endif  // FREHG2_BC_POLYGON_REGION_HPP
