#include "core/Simulation.hpp"

#include <string>
#include <type_traits>

int main()
{
    static_assert(std::is_same_v<frehg2::real, double>);
    static_assert(std::is_unsigned_v<frehg2::index_t>);

    frehg2::GridSpec grid;
    grid.nx = 4;
    grid.ny = 3;
    grid.nz = 2;
    grid.dx = 10.0;

    if (grid.nx != 4 || grid.ny != 3 || grid.nz != 2 || grid.dx != 10.0) {
        return 1;
    }

    frehg2::SimulationConfig config;
    config.input_path = "case.yaml";
    config.mode = frehg2::SimMode::SW_ONLY;
    config.boundary_codes.surface_water = {1, 2, 3, 4};
    config.boundary_codes.groundwater = {1, 2, 3, 4, 5, 6};

    if (std::string(frehg2::simModeName(config.mode)) != "surface_water") {
        return 1;
    }
    if (config.boundary_codes.surface_water.size() != 4) {
        return 1;
    }
    if (config.boundary_codes.groundwater.size() != 6) {
        return 1;
    }

    const frehg2::Simulation simulation(config);
    if (simulation.config().input_path != "case.yaml") {
        return 1;
    }
    if (simulation.run() != 0) {
        return 1;
    }

    const auto build_info = frehg2::getBuildInfo();
    (void)build_info;

    return 0;
}
