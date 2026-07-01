// Named performance regions for P21 timing breakdown (Task 21.3.1).
#ifndef FREHG2_PERF_PERF_REGIONS_HPP
#define FREHG2_PERF_PERF_REGIONS_HPP

namespace frehg2 {
namespace perf {

// Stable region ids written to simulation_summary.txt as perf_<name>_seconds.
enum class Region : unsigned char {
  SwAssembly = 0,
  SwKsp,
  SwUpdate,       // halo, velocity, rain/evap, geometry (SW step minus assembly/KSP)
  GwAssembly,
  GwKsp,
  GwUpdate,
  Coupling,
  Solute,
  Io,
  Polygon,
  Count
};

inline const char* regionName(Region r) {
  switch (r) {
    case Region::SwAssembly: return "sw_assembly";
    case Region::SwKsp: return "sw_ksp";
    case Region::SwUpdate: return "sw_update";
    case Region::GwAssembly: return "gw_assembly";
    case Region::GwKsp: return "gw_ksp";
    case Region::GwUpdate: return "gw_update";
    case Region::Coupling: return "coupling";
    case Region::Solute: return "solute";
    case Region::Io: return "io";
    case Region::Polygon: return "polygon";
    default: return "unknown";
  }
}

}  // namespace perf
}  // namespace frehg2

#endif  // FREHG2_PERF_PERF_REGIONS_HPP
