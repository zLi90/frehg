#ifndef FREHG2_MONITORING_LINE_FLUX_HPP
#define FREHG2_MONITORING_LINE_FLUX_HPP

#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "monitoring/MonitorSpec.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

namespace frehg2 {

// Trapezoidal integration of (velocity · line_normal) * dl along a segment (2D surface or GW face).
double integrateLineFlux(const LineFluxSpec& line, const SweSolver* swe, const ReSolver* re,
                         const Grid& swe_grid, const Grid& gw_grid, int gnx, int gny, double dx,
                         double dy, double x0, double y0, const MpiComm* mc);

}  // namespace frehg2

#endif  // FREHG2_MONITORING_LINE_FLUX_HPP
