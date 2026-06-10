#ifndef FREHG2_CORE_SIMULATION_HPP
#define FREHG2_CORE_SIMULATION_HPP

#include "core/types.hpp"

#include <string>

namespace frehg2 {

class Simulation {
public:
    explicit Simulation(SimulationConfig config);

    int run() const;

    const SimulationConfig& config() const noexcept;

private:
    SimulationConfig config_;
};

SimulationConfig makeConfigFromCommandLine(int argc, char** argv);
std::string helpText(const char* executable_name);

}  // namespace frehg2

#endif  // FREHG2_CORE_SIMULATION_HPP
