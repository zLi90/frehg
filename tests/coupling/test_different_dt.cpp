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
    const frehg2::Coupling coupling;

    frehg2::AsyncCouplingClock async_clock;
    async_clock.async = true;
    const auto async_steps = coupling.groundwaterStepsToSurfaceTime(async_clock, 10.0, 2.5);
    if (async_steps != 4 || !near(async_clock.surface_time, 10.0, 1.0e-12) ||
        !near(async_clock.groundwater_time, 10.0, 1.0e-12)) {
        std::cerr << "async coupling did not take the expected groundwater substeps\n";
        return 1;
    }

    const auto lagged_steps = coupling.groundwaterStepsToSurfaceTime(async_clock, 3.0, 2.5);
    if (lagged_steps != 1 || !near(async_clock.surface_time, 13.0, 1.0e-12) ||
        !near(async_clock.groundwater_time, 12.5, 1.0e-12)) {
        std::cerr << "async groundwater clock should remain within one GW dt of surface time\n";
        return 1;
    }

    frehg2::AsyncCouplingClock sync_clock;
    sync_clock.async = false;
    const auto sync_steps = coupling.groundwaterStepsToSurfaceTime(sync_clock, 10.0, 2.5);
    if (sync_steps != 1 || !near(sync_clock.groundwater_time, sync_clock.surface_time, 1.0e-12)) {
        std::cerr << "sync coupling should pin groundwater time to surface time\n";
        return 1;
    }

    return 0;
}
