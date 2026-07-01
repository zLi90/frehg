// Monitor probe / line-flux specifications (P15.3.1 / 15.3.2).
#ifndef FREHG2_MONITORING_MONITOR_SPEC_HPP
#define FREHG2_MONITORING_MONITOR_SPEC_HPP

#include <array>
#include <string>
#include <vector>

namespace frehg2 {

struct ProbeSpec {
  std::string name;
  // Physical coordinates (m); z used for subsurface probes.
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  // Resolved global cell indices after build().
  int gi = 0;
  int gj = 0;
  int gk = 0;
  int owner_rank = 0;
  int local_i = 0;
  int local_j = 0;
  int local_k = 0;
  bool subsurface = false;
  bool indices_from_grid = false;
  bool coords_from_xyz = false;
  std::vector<std::string> fields;
};

struct LineFluxSpec {
  std::string name;
  std::array<double, 3> p0{{0.0, 0.0, 0.0}};
  std::array<double, 3> p1{{0.0, 0.0, 0.0}};
  // Velocity component for surface lines ("u","v") or Darcy flux ("qx","qy") for GW sections.
  std::string field = "u";
  bool subsurface = false;
};

struct MonitorBundle {
  std::vector<ProbeSpec> probes;
  std::vector<LineFluxSpec> lines;
};

}  // namespace frehg2

#endif  // FREHG2_MONITORING_MONITOR_SPEC_HPP
