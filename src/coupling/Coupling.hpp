#ifndef FREHG2_COUPLING_COUPLING_HPP
#define FREHG2_COUPLING_COUPLING_HPP

#include "core/define.hpp"

#include <cstddef>
#include <vector>

namespace frehg2 {

struct CouplingParameters {
    real min_depth = 1.0e-8;
    real water_content_saturation = 1.0;
};

struct CoupledColumn {
    real surface_depth = 0.0;
    real groundwater_head = 0.0;
    real dz_top = 1.0;
    real saturated_conductivity = 0.0;
    real face_area = 1.0;
};

struct AsyncCoupledColumn {
    real surface_depth = 0.0;
    real surface_eta = 0.0;
    real groundwater_storage_depth = 0.0;
    real groundwater_head = 0.0;
    real dz_top = 1.0;
    real saturated_conductivity = 0.0;
    real face_area = 1.0;
};

struct AsyncCouplingClock {
    bool async = true;
    real surface_time = 0.0;
    real groundwater_time = 0.0;
};

struct AsyncCouplingResult {
    std::vector<real> exchange_rate;
    real surface_volume_change = 0.0;
    real groundwater_volume_change = 0.0;
    real surface_time = 0.0;
    real groundwater_time = 0.0;
    std::size_t groundwater_steps = 0;
};

class Coupling {
public:
    explicit Coupling(CouplingParameters parameters = {});

    const CouplingParameters& parameters() const noexcept;

    real computeExchangeFlux(const CoupledColumn& column) const;
    real computeExchangeRate(const CoupledColumn& column) const;
    real limitInfiltrationBySurfaceWater(real exchange_rate, real surface_depth, real dt) const;
    real synchronizedTimeStep(real surface_dt, real groundwater_dt) const;
    std::size_t groundwaterStepsToSurfaceTime(
        AsyncCouplingClock& clock,
        real surface_dt,
        real groundwater_dt) const;

    void applyExchangeToSurface(
        std::vector<real>& surface_depth,
        std::vector<real>& surface_eta,
        const std::vector<real>& exchange_rate,
        real dt) const;

    AsyncCouplingResult asyncCoupling(
        std::vector<AsyncCoupledColumn>& columns,
        AsyncCouplingClock& clock,
        real surface_dt,
        real groundwater_dt,
        real rain_rate,
        real evap_rate = 0.0) const;

private:
    CouplingParameters parameters_;
};

}  // namespace frehg2

#endif  // FREHG2_COUPLING_COUPLING_HPP
