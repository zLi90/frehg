// Initial-condition field specifications (P14.3.1).
#ifndef FREHG2_IC_IC_SPEC_HPP
#define FREHG2_IC_IC_SPEC_HPP

#include <string>
#include <unordered_map>
#include <vector>

#include "frehg2/core/define.hpp"
#include "bc/Polygon.hpp"

namespace frehg2 {

enum class ICKind { Constant, Raster, Polygon, Formula, Restart };

// Which solver field an IC block targets.
enum class ICField {
  SurfaceEta,
  SurfaceU,
  SurfaceV,
  GroundwaterWc,
  GroundwaterHead,
  SoluteSurface,
  SoluteSubsurface
};

struct ICFieldSpec {
  ICKind kind = ICKind::Constant;
  real constant_value = 0.0;
  std::string file;
  std::string format = "auto";  // auto | list | raster | hdf5
  std::string dataset = "data"; // HDF5 dataset path
  std::string formula;
  real polygon_default = 0.0;
  // Polygon name -> prescribed value (uses initial_conditions.regions polygons).
  std::unordered_map<std::string, real> polygon_values;
  std::string restart_file;
  real restart_time = 0.0;
};

struct ICRegion {
  Polygon polygon;
  real value = 0.0;
};

// Parsed flexible-IC bundle for one simulation.
struct InitialConditionsConfig {
  bool use_restart = false;
  std::string restart_file;
  real restart_time = 0.0;

  ICFieldSpec surface_eta;
  ICFieldSpec surface_u;
  ICFieldSpec surface_v;
  ICFieldSpec groundwater_wc;
  ICFieldSpec groundwater_head;
  ICFieldSpec solute_surface;
  ICFieldSpec solute_subsurface;

  std::vector<ICRegion> regions;
};

}  // namespace frehg2

#endif  // FREHG2_IC_IC_SPEC_HPP
