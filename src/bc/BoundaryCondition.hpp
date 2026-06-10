#ifndef FREHG2_BC_BOUNDARY_CONDITION_HPP
#define FREHG2_BC_BOUNDARY_CONDITION_HPP

#include "bc/Polygon.hpp"

#include <vector>

namespace frehg2 {

enum class BoundaryConditionType {
    FIXED_WATER_LEVEL,
    FIXED_FLOW_RATE,
    FREE_OUTFLOW,
    ZERO_GRADIENT,
    TIDAL_WATER_LEVEL,
    GROUNDWATER_HEAD,
    ROOT_WATER_UPTAKE,
    GRAVITY_DRAINAGE
};

enum class SourceSinkType {
    RAINFALL,
    EVAPOTRANSPIRATION,
    INFLOW,
    PUMPING
};

struct PolygonBoundaryCondition {
    Polygon polygon;
    BoundaryConditionType type = BoundaryConditionType::ZERO_GRADIENT;
    real value = 0.0;
    real normal_x = 1.0;
    real normal_y = 0.0;
};

struct PolygonSourceSink {
    Polygon polygon;
    SourceSinkType type = SourceSinkType::RAINFALL;
    real rate = 0.0;
};

class BoundaryCondition {
public:
    static std::vector<int> markCells(
        const Grid& grid,
        const std::vector<PolygonBoundaryCondition>& conditions);

    static void applySurface(
        const Grid& grid,
        const std::vector<PolygonBoundaryCondition>& conditions,
        std::vector<real>& eta,
        const std::vector<real>& bed);

    static void applyGroundwaterHead(
        const Grid& grid,
        const std::vector<PolygonBoundaryCondition>& conditions,
        std::vector<real>& head);

    static bool applyGroundwaterTopFlux(
        const Grid& grid,
        const std::vector<PolygonBoundaryCondition>& conditions,
        std::vector<real>& top_flux);

    static void applySourceSink(
        const Grid& grid,
        const std::vector<PolygonSourceSink>& sources,
        std::vector<real>& values,
        real dt);
};

}  // namespace frehg2

#endif  // FREHG2_BC_BOUNDARY_CONDITION_HPP
