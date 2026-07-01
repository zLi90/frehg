#ifndef FREHG2_MONITORING_PROBE_LOCATOR_HPP
#define FREHG2_MONITORING_PROBE_LOCATOR_HPP

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "monitoring/MonitorSpec.hpp"

namespace frehg2 {

// Resolve probe global/local indices and owner rank from physical xyz or preset gi/gj/gk.
void buildProbeLocations(MonitorBundle& bundle, int gnx, int gny, int gnz, double dx, double dy,
                         double dz, double x0, double y0, double botz, const MpiComm* mc);

}  // namespace frehg2

#endif  // FREHG2_MONITORING_PROBE_LOCATOR_HPP
