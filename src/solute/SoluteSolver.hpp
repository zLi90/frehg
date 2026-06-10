#ifndef FREHG2_SOLUTE_SOLVER_HPP
#define FREHG2_SOLUTE_SOLVER_HPP

#include "core/define.hpp"

#include <vector>

namespace frehg2 {

struct SoluteParameters {
    real dt = 1.0;
    real dx = 1.0;
    real diffusion = 0.0;
    real min_concentration = 0.0;
    real max_concentration = 200.0;
    bool use_superbee = false;
};

struct SoluteExchange {
    real exchange_rate = 0.0;
    real interface_diffusion = 0.0;
    real distance = 1.0;
    real area = 1.0;
};

class SoluteSolver {
public:
    explicit SoluteSolver(SoluteParameters parameters);

    const SoluteParameters& parameters() const noexcept;

    static real superbee(real sp, real sc, real sm, real velocity, real delta, real dt);
    static real upwindFaceValue(real left, real right, real face_flux);
    static real exchangeFlux(real surface_concentration, real subsurface_concentration,
                             const SoluteExchange& exchange);

    std::vector<real> advect1D(
        const std::vector<real>& concentration,
        const std::vector<real>& face_flux,
        const std::vector<real>& volume) const;

    std::vector<real> diffuse1D(
        const std::vector<real>& concentration,
        const std::vector<real>& volume) const;

    void applyConservativeExchange(
        real& surface_mass,
        real& subsurface_mass,
        real surface_concentration,
        real subsurface_concentration,
        const SoluteExchange& exchange) const;

    void limit(std::vector<real>& concentration) const;

private:
    SoluteParameters parameters_;
};

}  // namespace frehg2

#endif  // FREHG2_SOLUTE_SOLVER_HPP
