#include "solute/SoluteSolver.hpp"

#include <cmath>
#include <iostream>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

}  // namespace

int main()
{
    frehg2::SoluteParameters parameters;
    parameters.dt = 2.0;
    frehg2::SoluteSolver solver(parameters);

    frehg2::SoluteExchange seepage;
    seepage.exchange_rate = 0.01;
    seepage.interface_diffusion = 0.0;
    seepage.distance = 0.5;
    seepage.area = 10.0;

    double surface_mass = 5.0;
    double subsurface_mass = 20.0;
    const double initial_total = surface_mass + subsurface_mass;
    solver.applyConservativeExchange(surface_mass, subsurface_mass, 1.0, 4.0, seepage);

    const double expected_transfer = seepage.exchange_rate * 4.0 * seepage.area * parameters.dt;
    if (!near(surface_mass, 5.0 + expected_transfer, 1.0e-12) ||
        !near(subsurface_mass, 20.0 - expected_transfer, 1.0e-12)) {
        std::cerr << "seepage scalar exchange used the wrong concentration\n";
        return 1;
    }
    if (!near(initial_total, surface_mass + subsurface_mass, 1.0e-12)) {
        std::cerr << "solute exchange did not conserve total mass\n";
        return 1;
    }

    frehg2::SoluteExchange infiltration;
    infiltration.exchange_rate = -0.01;
    infiltration.interface_diffusion = 0.0;
    infiltration.distance = 0.5;
    infiltration.area = 10.0;
    const auto flux = frehg2::SoluteSolver::exchangeFlux(1.0, 4.0, infiltration);
    if (!near(flux, -0.1, 1.0e-12)) {
        std::cerr << "infiltration scalar flux should use surface concentration\n";
        return 1;
    }

    frehg2::SoluteExchange diffusion;
    diffusion.exchange_rate = 0.0;
    diffusion.interface_diffusion = 0.2;
    diffusion.distance = 0.5;
    diffusion.area = 10.0;
    if (!near(frehg2::SoluteSolver::exchangeFlux(1.0, 4.0, diffusion), 24.0, 1.0e-12)) {
        std::cerr << "diffusive scalar exchange mismatch\n";
        return 1;
    }
    return 0;
}
