// YAML parsing for flexible initial conditions (P14.3.2 / 14.3.4).
#ifndef FREHG2_IC_IC_CONFIG_HPP
#define FREHG2_IC_IC_CONFIG_HPP

#include "ic/ICSpec.hpp"

namespace frehg2 {

class Config;

// default_gw_wc is used when groundwater.wc is absent (typically soil.theta_r).
InitialConditionsConfig parseInitialConditions(const Config& config, real default_gw_wc = 0.0);

}  // namespace frehg2

#endif  // FREHG2_IC_IC_CONFIG_HPP
