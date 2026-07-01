// YAML parsing for the non-uniform soil schema (P13.3.2).
//
// `soil.types` is an ordered list of soil CLASSES (class id == list position). Each class carries
// the van Genuchten / Ksat tuple; `specific_storage`, `use_vg`, `use_mvg`, and `air_entry_value`
// default to the global `groundwater.*` settings (so class 0 reproduces the uniform P5 rp.soil
// exactly) but may be overridden per class.
#ifndef FREHG2_SOIL_SOIL_CONFIG_HPP
#define FREHG2_SOIL_SOIL_CONFIG_HPP

#include <vector>

#include "frehg2/re/SoilParams.hpp"

namespace frehg2 {

class Config;

// Parse all `soil.types[]` classes. Throws if `soil.types` is missing/empty or any class omits a
// required field. Returns at least one class.
std::vector<SoilParams> parseSoilClasses(const Config& config);

}  // namespace frehg2

#endif  // FREHG2_SOIL_SOIL_CONFIG_HPP
