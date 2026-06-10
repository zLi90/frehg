#include "coupling/Coupling.hpp"

#include <algorithm>
#include <stdexcept>

namespace frehg2 {

Coupling::Coupling(CouplingParameters parameters)
    : parameters_(parameters)
{
    if (parameters_.water_content_saturation <= 0.0) {
        throw std::invalid_argument("saturated water content must be positive");
    }
}

const CouplingParameters& Coupling::parameters() const noexcept
{
    return parameters_;
}

real Coupling::computeExchangeFlux(const CoupledColumn& column) const
{
    if (column.dz_top <= 0.0 || column.face_area <= 0.0) {
        throw std::invalid_argument("coupled column geometry must be positive");
    }
    const real delta = 0.5 * column.dz_top;
    const real rate =
        column.saturated_conductivity *
        ((column.groundwater_head - column.surface_depth) / delta - 1.0);
    return column.face_area * rate;
}

real Coupling::computeExchangeRate(const CoupledColumn& column) const
{
    return computeExchangeFlux(column) / column.face_area;
}

real Coupling::limitInfiltrationBySurfaceWater(real exchange_rate, real surface_depth, real dt) const
{
    if (dt <= 0.0) {
        throw std::invalid_argument("coupling timestep must be positive");
    }
    if (exchange_rate >= 0.0) {
        return exchange_rate;
    }
    if (surface_depth <= parameters_.min_depth) {
        return 0.0;
    }
    const real max_infiltration = surface_depth / (dt * parameters_.water_content_saturation);
    return std::max(exchange_rate, -max_infiltration);
}

real Coupling::synchronizedTimeStep(real surface_dt, real groundwater_dt) const
{
    if (surface_dt <= 0.0 || groundwater_dt <= 0.0) {
        throw std::invalid_argument("solver timesteps must be positive");
    }
    return std::min(surface_dt, groundwater_dt);
}

std::size_t Coupling::groundwaterStepsToSurfaceTime(
    AsyncCouplingClock& clock,
    real surface_dt,
    real groundwater_dt) const
{
    if (surface_dt <= 0.0 || groundwater_dt <= 0.0) {
        throw std::invalid_argument("solver timesteps must be positive");
    }

    clock.surface_time += surface_dt;
    if (!clock.async) {
        clock.groundwater_time = clock.surface_time;
        return 1;
    }

    std::size_t steps = 0;
    constexpr real TIME_EPS = 1.0e-12;
    while (clock.groundwater_time + groundwater_dt <= clock.surface_time + TIME_EPS) {
        clock.groundwater_time += groundwater_dt;
        ++steps;
    }
    return steps;
}

void Coupling::applyExchangeToSurface(
    std::vector<real>& surface_depth,
    std::vector<real>& surface_eta,
    const std::vector<real>& exchange_rate,
    real dt) const
{
    if (dt <= 0.0) {
        throw std::invalid_argument("coupling timestep must be positive");
    }
    if (surface_depth.size() != surface_eta.size() || surface_depth.size() != exchange_rate.size()) {
        throw std::invalid_argument("surface and exchange arrays must have matching sizes");
    }

    for (std::size_t i = 0; i < surface_depth.size(); ++i) {
        const real limited_rate =
            limitInfiltrationBySurfaceWater(exchange_rate[i], surface_depth[i], dt);
        const real delta_depth = limited_rate * dt * parameters_.water_content_saturation;
        surface_depth[i] = std::max<real>(0.0, surface_depth[i] + delta_depth);
        surface_eta[i] += delta_depth;
    }
}

AsyncCouplingResult Coupling::asyncCoupling(
    std::vector<AsyncCoupledColumn>& columns,
    AsyncCouplingClock& clock,
    real surface_dt,
    real groundwater_dt,
    real rain_rate,
    real evap_rate) const
{
    if (surface_dt <= 0.0 || groundwater_dt <= 0.0) {
        throw std::invalid_argument("solver timesteps must be positive");
    }
    AsyncCouplingResult result;
    result.exchange_rate.assign(columns.size(), 0.0);

    const real requested_rain_delta = (rain_rate - evap_rate) * surface_dt;
    for (auto& column : columns) {
        if (column.dz_top <= 0.0 || column.face_area <= 0.0) {
            throw std::invalid_argument("coupled column geometry must be positive");
        }
        const real rain_delta =
            requested_rain_delta < 0.0
                ? std::max(requested_rain_delta, -column.surface_depth)
                : requested_rain_delta;
        column.surface_depth += rain_delta;
        column.surface_eta += rain_delta;
        result.surface_volume_change += rain_delta * column.face_area;
    }

    for (std::size_t i = 0; i < columns.size(); ++i) {
        auto& column = columns[i];
        const CoupledColumn exchange_column{
            column.surface_depth,
            column.groundwater_head,
            column.dz_top,
            column.saturated_conductivity,
            column.face_area,
        };
        const real raw_rate = computeExchangeRate(exchange_column);
        const real limited_rate =
            limitInfiltrationBySurfaceWater(raw_rate, column.surface_depth, surface_dt);
        const real delta_depth = limited_rate * surface_dt * parameters_.water_content_saturation;
        const real exchange_volume = delta_depth * column.face_area;

        column.surface_depth = std::max<real>(0.0, column.surface_depth + delta_depth);
        column.surface_eta += delta_depth;
        column.groundwater_storage_depth -= delta_depth;

        result.exchange_rate[i] = limited_rate;
        result.surface_volume_change += exchange_volume;
        result.groundwater_volume_change -= exchange_volume;
    }

    result.groundwater_steps = groundwaterStepsToSurfaceTime(clock, surface_dt, groundwater_dt);
    result.surface_time = clock.surface_time;
    result.groundwater_time = clock.groundwater_time;
    return result;
}

}  // namespace frehg2
