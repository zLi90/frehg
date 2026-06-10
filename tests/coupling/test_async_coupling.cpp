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
    parameters.water_content_saturation = 1.0;
    const frehg2::Coupling coupling(parameters);

    std::vector<frehg2::AsyncCoupledColumn> columns(1);
    columns[0].surface_depth = 0.0;
    columns[0].surface_eta = 0.0;
    columns[0].groundwater_storage_depth = 0.0;
    columns[0].groundwater_head = -0.05;
    columns[0].dz_top = 0.10;
    columns[0].saturated_conductivity = 1.0e-4;
    columns[0].face_area = 2.0;

    frehg2::AsyncCouplingClock clock;
    clock.async = true;

    const auto result = coupling.asyncCoupling(
        columns,
        clock,
        10.0,
        5.0,
        0.001);

    if (!(result.exchange_rate[0] < 0.0)) {
        std::cerr << "rain-created ponding should infiltrate into groundwater\n";
        return 1;
    }
    if (columns[0].surface_depth < -1.0e-14) {
        std::cerr << "async coupling drained more surface water than was available\n";
        return 1;
    }
    if (!(columns[0].groundwater_storage_depth > 0.0)) {
        std::cerr << "groundwater did not gain infiltrated water\n";
        return 1;
    }
    if (!near(result.surface_volume_change + result.groundwater_volume_change,
              0.001 * 10.0 * columns[0].face_area,
              1.0e-12)) {
        std::cerr << "rain to ponding to infiltration mass balance failed\n";
        return 1;
    }
    if (result.groundwater_steps != 2 || !near(result.groundwater_time, 10.0, 1.0e-12)) {
        std::cerr << "groundwater clock did not advance to the surface step\n";
        return 1;
    }
    return 0;
}
