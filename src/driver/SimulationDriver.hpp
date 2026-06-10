#ifndef FREHG2_DRIVER_SIMULATION_DRIVER_HPP
#define FREHG2_DRIVER_SIMULATION_DRIVER_HPP

#include <filesystem>

namespace frehg2 {

class SimulationDriver {
public:
    explicit SimulationDriver(std::filesystem::path config_path);

    int run() const;

private:
    std::filesystem::path config_path_;

    int runBenchmarkSurfaceWater() const;
    int runBenchmarkGroundwater() const;
    int runCoupledSmoke() const;
};

}  // namespace frehg2

#endif  // FREHG2_DRIVER_SIMULATION_DRIVER_HPP
