#ifndef FREHG2_MONITORING_MONITOR_CONFIG_HPP
#define FREHG2_MONITORING_MONITOR_CONFIG_HPP

#include "monitoring/MonitorSpec.hpp"

namespace frehg2 {

class Config;

MonitorBundle parseMonitors(const Config& config, bool sw_enabled, bool gw_enabled);

}  // namespace frehg2

#endif  // FREHG2_MONITORING_MONITOR_CONFIG_HPP
