#include "bc/PolygonRegion.hpp"

#include <stdexcept>

namespace frehg2 {

BcKind bcKindFromString(const std::string& s) {
  if (s == "bc_discharge") return BcKind::Discharge;
  if (s == "bc_depth") return BcKind::Depth;
  if (s == "bc_critical") return BcKind::Critical;
  throw std::runtime_error("PolygonRegion: unknown boundary type '" + s +
                           "' (expected bc_discharge | bc_depth | bc_critical)");
}

SourceKind sourceKindFromString(const std::string& s) {
  if (s == "inflow_rate") return SourceKind::InflowRate;
  if (s == "rainfall_rate") return SourceKind::RainfallRate;
  if (s == "extraction_well") return SourceKind::ExtractionWell;
  throw std::runtime_error("PolygonRegion: unknown source type '" + s +
                           "' (expected inflow_rate | rainfall_rate | extraction_well)");
}

const char* toString(BcKind k) {
  switch (k) {
    case BcKind::Discharge: return "bc_discharge";
    case BcKind::Depth: return "bc_depth";
    case BcKind::Critical: return "bc_critical";
  }
  return "bc_discharge";
}

const char* toString(SourceKind k) {
  switch (k) {
    case SourceKind::InflowRate: return "inflow_rate";
    case SourceKind::RainfallRate: return "rainfall_rate";
    case SourceKind::ExtractionWell: return "extraction_well";
  }
  return "inflow_rate";
}

bool isBcTypeString(const std::string& s) {
  return s == "bc_discharge" || s == "bc_depth" || s == "bc_critical";
}

bool isSourceTypeString(const std::string& s) {
  return s == "inflow_rate" || s == "rainfall_rate" || s == "extraction_well";
}

}  // namespace frehg2
