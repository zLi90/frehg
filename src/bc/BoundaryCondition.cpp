#include "bc/BoundaryCondition.hpp"

#include <algorithm>
#include <stdexcept>

namespace frehg2 {

namespace {

int code(BoundaryConditionType type)
{
    switch (type) {
    case BoundaryConditionType::FIXED_WATER_LEVEL:
        return 1;
    case BoundaryConditionType::FIXED_FLOW_RATE:
        return 2;
    case BoundaryConditionType::FREE_OUTFLOW:
        return 3;
    case BoundaryConditionType::ZERO_GRADIENT:
        return 4;
    case BoundaryConditionType::TIDAL_WATER_LEVEL:
        return 5;
    case BoundaryConditionType::GROUNDWATER_HEAD:
        return 6;
    case BoundaryConditionType::ROOT_WATER_UPTAKE:
        return 7;
    case BoundaryConditionType::GRAVITY_DRAINAGE:
        return 8;
    }
    return 0;
}

real sourceSign(SourceSinkType type)
{
    switch (type) {
    case SourceSinkType::RAINFALL:
    case SourceSinkType::INFLOW:
        return 1.0;
    case SourceSinkType::EVAPOTRANSPIRATION:
    case SourceSinkType::PUMPING:
        return -1.0;
    }
    return 0.0;
}

}  // namespace

std::vector<int> BoundaryCondition::markCells(
    const Grid& grid,
    const std::vector<PolygonBoundaryCondition>& conditions)
{
    std::vector<int> markers(static_cast<std::size_t>(grid.nSurfaceCell()), 0);
    for (const auto& condition : conditions) {
        for (const auto cell : condition.polygon.selectSurfaceCells(grid)) {
            markers[cell] = code(condition.type);
        }
    }
    return markers;
}

void BoundaryCondition::applySurface(
    const Grid& grid,
    const std::vector<PolygonBoundaryCondition>& conditions,
    std::vector<real>& eta,
    const std::vector<real>& bed)
{
    if (eta.size() < static_cast<std::size_t>(grid.nSurfaceCell()) ||
        bed.size() < static_cast<std::size_t>(grid.nSurfaceCell())) {
        throw std::invalid_argument("surface arrays do not cover the grid");
    }

    for (const auto& condition : conditions) {
        const auto cells = condition.polygon.selectSurfaceCells(grid);
        for (const auto cell : cells) {
            switch (condition.type) {
            case BoundaryConditionType::FIXED_WATER_LEVEL:
            case BoundaryConditionType::TIDAL_WATER_LEVEL:
                eta[cell] = std::max(bed[cell], condition.value);
                break;
            case BoundaryConditionType::FIXED_FLOW_RATE:
            case BoundaryConditionType::FREE_OUTFLOW:
            case BoundaryConditionType::ZERO_GRADIENT:
            case BoundaryConditionType::GROUNDWATER_HEAD:
            case BoundaryConditionType::ROOT_WATER_UPTAKE:
            case BoundaryConditionType::GRAVITY_DRAINAGE:
                break;
            }
        }
    }
}

void BoundaryCondition::applyGroundwaterHead(
    const Grid& grid,
    const std::vector<PolygonBoundaryCondition>& conditions,
    std::vector<real>& head)
{
    if (head.size() < static_cast<std::size_t>(grid.nCell())) {
        throw std::invalid_argument("groundwater head array does not cover the grid");
    }

    for (const auto& condition : conditions) {
        if (condition.type != BoundaryConditionType::GROUNDWATER_HEAD) {
            continue;
        }
        for (const auto cell : condition.polygon.selectGroundwaterCells(grid)) {
            head[cell] = condition.value;
        }
    }
}

bool BoundaryCondition::applyGroundwaterTopFlux(
    const Grid& grid,
    const std::vector<PolygonBoundaryCondition>& conditions,
    std::vector<real>& top_flux)
{
    if (top_flux.size() < static_cast<std::size_t>(grid.nSurfaceCell())) {
        throw std::invalid_argument("groundwater top-flux array does not cover the surface grid");
    }

    bool applied = false;
    for (const auto& condition : conditions) {
        if (condition.type != BoundaryConditionType::FIXED_FLOW_RATE) {
            continue;
        }
        for (const auto cell : condition.polygon.selectSurfaceCells(grid)) {
            top_flux[cell] = condition.value;
            applied = true;
        }
    }
    return applied;
}

void BoundaryCondition::applySourceSink(
    const Grid& grid,
    const std::vector<PolygonSourceSink>& sources,
    std::vector<real>& values,
    real dt)
{
    if (dt <= 0.0) {
        throw std::invalid_argument("source/sink timestep must be positive");
    }
    if (values.size() < static_cast<std::size_t>(grid.nSurfaceCell())) {
        throw std::invalid_argument("source/sink target array does not cover the surface grid");
    }

    for (const auto& source : sources) {
        const real delta = sourceSign(source.type) * source.rate * dt;
        for (const auto cell : source.polygon.selectSurfaceCells(grid)) {
            values[cell] = std::max<real>(0.0, values[cell] + delta);
        }
    }
}

}  // namespace frehg2
