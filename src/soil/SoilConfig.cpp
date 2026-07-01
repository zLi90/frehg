#include "soil/SoilConfig.hpp"

#include <stdexcept>

#include "io/Config.hpp"

namespace frehg2 {

namespace {
SoilParams parseOne(const YAML::Node& st, const Config& config) {
  if (!st["theta_s"] || !st["theta_r"] || !st["vg"] || !st["k_sat"])
    throw std::runtime_error(
        "parseSoilClasses: a soil.types entry is missing theta_s/theta_r/vg/k_sat");
  SoilParams sp;
  sp.theta_s = st["theta_s"].as<double>();
  sp.theta_r = st["theta_r"].as<double>();
  sp.alpha = st["vg"]["alpha"].as<double>();
  sp.n = st["vg"]["n"].as<double>();
  sp.Ks_x = st["k_sat"]["x"].as<double>();
  sp.Ks_y = st["k_sat"]["y"].as<double>();
  sp.Ks_z = st["k_sat"]["z"].as<double>();
  sp.Ss = config.getOr<double>("groundwater.specific_storage", 1.0e-5);
  sp.use_vg = config.getOr<bool>("groundwater.use_vg", true);
  sp.use_mvg = config.getOr<bool>("groundwater.use_mvg", false);
  sp.air_entry = config.getOr<double>("groundwater.air_entry_value", -0.02);
  if (st["specific_storage"]) sp.Ss = st["specific_storage"].as<double>();
  if (st["use_vg"]) sp.use_vg = st["use_vg"].as<bool>();
  if (st["use_mvg"]) sp.use_mvg = st["use_mvg"].as<bool>();
  if (st["air_entry_value"]) sp.air_entry = st["air_entry_value"].as<double>();
  return sp;
}
}  // namespace

std::vector<SoilParams> parseSoilClasses(const Config& config) {
  const YAML::Node soil_types = config.root()["soil"]["types"];
  if (!soil_types || !soil_types.IsSequence() || soil_types.size() == 0)
    throw std::runtime_error("parseSoilClasses: soil.types is missing or empty");
  std::vector<SoilParams> classes;
  classes.reserve(soil_types.size());
  for (size_t c = 0; c < soil_types.size(); ++c)
    classes.push_back(parseOne(soil_types[c], config));
  return classes;
}

}  // namespace frehg2
