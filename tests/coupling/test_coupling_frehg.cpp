#include "coupling/Coupling.hpp"

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
    frehg2::Coupling coupling;
    const frehg2::CoupledColumn column{
        0.10,      // surface depth
        -0.05,     // top-cell groundwater pressure head
        0.01,      // top layer thickness
        2.89e-6,   // saturated Kz
        4.0,       // face area
    };

    const auto rate = coupling.computeExchangeRate(column);
    const auto expected_rate = 2.89e-6 * (((-0.05) - 0.10) / (0.5 * 0.01) - 1.0);
    if (!near(rate, expected_rate, 1.0e-18)) {
        std::cerr << "Frehg exchange rate mismatch\n";
        return 1;
    }
    if (!(rate < 0.0)) {
        std::cerr << "ponded surface water should infiltrate with negative exchange\n";
        return 1;
    }

    const auto flux = coupling.computeExchangeFlux(column);
    if (!near(flux, expected_rate * column.face_area, 1.0e-18)) {
        std::cerr << "Frehg exchange flux mismatch\n";
        return 1;
    }
    return 0;
}
