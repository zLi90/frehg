#include "coupling/Coupling.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

}  // namespace

int main()
{
    frehg2::CouplingParameters parameters;
    parameters.min_depth = 1.0e-8;
    parameters.water_content_saturation = 0.33;
    frehg2::Coupling coupling(parameters);

    const auto dt = coupling.synchronizedTimeStep(2.0, 5.0);
    if (!near(dt, 2.0, 0.0)) {
        std::cerr << "sync timestep should use the smaller solver dt\n";
        return 1;
    }

    std::vector<frehg2::real> depth{0.02};
    std::vector<frehg2::real> eta{1.02};
    const std::vector<frehg2::real> exchange_rate{-0.10};

    const auto initial_surface = depth[0];
    const auto limited_rate = coupling.limitInfiltrationBySurfaceWater(exchange_rate[0], depth[0], dt);
    coupling.applyExchangeToSurface(depth, eta, exchange_rate, dt);
    const auto gained_groundwater = -limited_rate * dt * parameters.water_content_saturation;

    if (depth[0] < -1.0e-14) {
        std::cerr << "surface depth became negative\n";
        return 1;
    }
    if (!near(initial_surface, depth[0] + gained_groundwater, 1.0e-12)) {
        std::cerr << "coupled surface/subsurface mass balance failed\n";
        return 1;
    }
    if (!near(eta[0], 1.02 + (depth[0] - initial_surface), 1.0e-12)) {
        std::cerr << "surface elevation was not adjusted with depth\n";
        return 1;
    }
    return 0;
}
