#include "solute/SoluteSolver.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

double totalMass(const std::vector<frehg2::real>& concentration, const std::vector<frehg2::real>& volume)
{
    double total = 0.0;
    for (std::size_t i = 0; i < concentration.size(); ++i) {
        total += concentration[i] * volume[i];
    }
    return total;
}

}  // namespace

int main()
{
    frehg2::SoluteParameters parameters;
    parameters.dt = 0.1;
    parameters.dx = 1.0;
    parameters.diffusion = 0.1;
    frehg2::SoluteSolver solver(parameters);

    const std::vector<frehg2::real> volume{1.0, 1.0, 1.0, 1.0, 1.0};
    const std::vector<frehg2::real> initial{0.0, 0.0, 10.0, 0.0, 0.0};
    const std::vector<frehg2::real> face_flux{0.0, 0.0, 0.4, 0.4, 0.0, 0.0};

    const auto advected = solver.advect1D(initial, face_flux, volume);
    if (!(advected[2] < initial[2] && advected[3] > initial[3])) {
        std::cerr << "upwind advection did not move the scalar pulse downstream\n";
        return 1;
    }
    if (!near(totalMass(initial, volume), totalMass(advected, volume), 1.0e-12)) {
        std::cerr << "solute advection did not conserve mass\n";
        return 1;
    }

    const auto diffused = solver.diffuse1D(initial, volume);
    if (!(diffused[1] > initial[1] && diffused[3] > initial[3] && diffused[2] < initial[2])) {
        std::cerr << "diffusion did not smooth the scalar pulse\n";
        return 1;
    }
    if (!near(totalMass(initial, volume), totalMass(diffused, volume), 1.0e-12)) {
        std::cerr << "solute diffusion did not conserve mass\n";
        return 1;
    }

    const auto tvd = frehg2::SoluteSolver::superbee(2.0, 1.0, 0.0, 0.25, 1.0, 0.1);
    if (!(tvd > 1.0 && tvd < 2.0)) {
        std::cerr << "superbee limiter returned an unexpected value\n";
        return 1;
    }
    return 0;
}
