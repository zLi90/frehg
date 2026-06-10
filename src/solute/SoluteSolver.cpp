#include "solute/SoluteSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace frehg2 {

SoluteSolver::SoluteSolver(SoluteParameters parameters)
    : parameters_(parameters)
{
    if (parameters_.dt <= 0.0 || parameters_.dx <= 0.0) {
        throw std::invalid_argument("solute timestep and grid spacing must be positive");
    }
    if (parameters_.min_concentration > parameters_.max_concentration) {
        throw std::invalid_argument("solute concentration limits are invalid");
    }
}

const SoluteParameters& SoluteSolver::parameters() const noexcept
{
    return parameters_;
}

real SoluteSolver::superbee(real sp, real sc, real sm, real velocity, real delta, real dt)
{
    if (delta <= 0.0 || dt <= 0.0) {
        throw std::invalid_argument("TVD spacing and timestep must be positive");
    }

    const real coef = std::abs(velocity) * dt / delta;
    real phi = 0.0;
    if (sp != sc) {
        const real r = (sc - sm) / (sp - sc);
        real r1 = 1.0;
        real r2 = 2.0;
        if (2.0 * r < 1.0) {
            r1 = 2.0 * r;
        }
        if (r < 2.0) {
            r2 = r;
        }
        if (r1 > 0.0 || r2 > 0.0) {
            phi = std::max(r1, r2);
        }
    }
    return sc + 0.5 * phi * (1.0 - coef) * (sp - sc);
}

real SoluteSolver::upwindFaceValue(real left, real right, real face_flux)
{
    return face_flux >= 0.0 ? left : right;
}

real SoluteSolver::exchangeFlux(
    real surface_concentration,
    real subsurface_concentration,
    const SoluteExchange& exchange)
{
    if (exchange.distance <= 0.0 || exchange.area <= 0.0) {
        throw std::invalid_argument("solute exchange geometry must be positive");
    }

    const real advective =
        exchange.exchange_rate > 0.0 ? exchange.exchange_rate * subsurface_concentration
                                     : exchange.exchange_rate * surface_concentration;
    const real diffusive =
        (2.0 * exchange.interface_diffusion / exchange.distance) *
        (subsurface_concentration - surface_concentration);
    return (advective + diffusive) * exchange.area;
}

std::vector<real> SoluteSolver::advect1D(
    const std::vector<real>& concentration,
    const std::vector<real>& face_flux,
    const std::vector<real>& volume) const
{
    const auto n = concentration.size();
    if (face_flux.size() != n + 1 || volume.size() != n) {
        throw std::invalid_argument("solute advection arrays have inconsistent sizes");
    }

    std::vector<real> mass(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (volume[i] <= 0.0) {
            throw std::invalid_argument("solute control-volume size must be positive");
        }
        mass[i] = concentration[i] * volume[i];
    }

    for (std::size_t face = 1; face < n; ++face) {
        real face_value =
            upwindFaceValue(concentration[face - 1], concentration[face], face_flux[face]);
        if (parameters_.use_superbee && face_flux[face] > 0.0) {
            const real left_left = face >= 2 ? concentration[face - 2] : concentration[face - 1];
            face_value = superbee(concentration[face], concentration[face - 1], left_left,
                                  face_flux[face], parameters_.dx, parameters_.dt);
        } else if (parameters_.use_superbee && face_flux[face] < 0.0) {
            const real right_right = face + 1 < n ? concentration[face + 1] : concentration[face];
            face_value = superbee(concentration[face - 1], concentration[face], right_right,
                                  face_flux[face], parameters_.dx, parameters_.dt);
        }
        const real transported_mass = face_flux[face] * face_value * parameters_.dt;
        mass[face - 1] -= transported_mass;
        mass[face] += transported_mass;
    }

    std::vector<real> updated(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        updated[i] = mass[i] / volume[i];
    }
    limit(updated);
    return updated;
}

std::vector<real> SoluteSolver::diffuse1D(
    const std::vector<real>& concentration,
    const std::vector<real>& volume) const
{
    const auto n = concentration.size();
    if (volume.size() != n) {
        throw std::invalid_argument("solute diffusion arrays have inconsistent sizes");
    }

    std::vector<real> mass(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (volume[i] <= 0.0) {
            throw std::invalid_argument("solute control-volume size must be positive");
        }
        mass[i] = concentration[i] * volume[i];
    }

    for (std::size_t face = 1; face < n; ++face) {
        const real flux =
            parameters_.diffusion * (concentration[face] - concentration[face - 1]) /
            parameters_.dx;
        const real transported_mass = flux * parameters_.dt;
        mass[face - 1] += transported_mass;
        mass[face] -= transported_mass;
    }

    std::vector<real> updated(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        updated[i] = mass[i] / volume[i];
    }
    limit(updated);
    return updated;
}

void SoluteSolver::applyConservativeExchange(
    real& surface_mass,
    real& subsurface_mass,
    real surface_concentration,
    real subsurface_concentration,
    const SoluteExchange& exchange) const
{
    const real mass_flux = exchangeFlux(surface_concentration, subsurface_concentration, exchange);
    const real transferred_mass = mass_flux * parameters_.dt;
    surface_mass += transferred_mass;
    subsurface_mass -= transferred_mass;
}

void SoluteSolver::limit(std::vector<real>& concentration) const
{
    for (auto& value : concentration) {
        value = std::max(parameters_.min_concentration,
                         std::min(value, parameters_.max_concentration));
    }
}

}  // namespace frehg2
