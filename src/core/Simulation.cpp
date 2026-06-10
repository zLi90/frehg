#include "core/Simulation.hpp"

#include <stdexcept>
#include <utility>

namespace frehg2 {

BuildInfo getBuildInfo()
{
    return BuildInfo{
#ifdef USE_KOKKOS
        true,
#else
        false,
#endif
#ifdef USE_MPI
        true,
#else
        false,
#endif
#ifdef USE_PETSC
        true,
#else
        false,
#endif
#ifdef USE_HDF5
        true,
#else
        false,
#endif
#ifdef USE_YAML_CPP
        true,
#else
        false,
#endif
    };
}

const char* simModeName(SimMode mode)
{
    switch (mode) {
    case SimMode::SW_ONLY:
        return "surface_water";
    case SimMode::GW_ONLY:
        return "groundwater";
    case SimMode::COUPLED:
        return "coupled";
    case SimMode::SOLUTE:
        return "solute";
    }

    return "unknown";
}

Simulation::Simulation(SimulationConfig config)
    : config_(std::move(config))
{
}

int Simulation::run() const
{
    return config_.input_path.empty() ? 1 : 0;
}

const SimulationConfig& Simulation::config() const noexcept
{
    return config_;
}

SimulationConfig makeConfigFromCommandLine(int argc, char** argv)
{
    if (argc != 2) {
        throw std::invalid_argument("expected exactly one configuration file path");
    }

    const std::string argument = argv[1] == nullptr ? std::string{} : std::string(argv[1]);
    if (argument == "--help" || argument == "-h") {
        throw std::invalid_argument("help requested");
    }

    SimulationConfig config;
    config.input_path = argument;
    return config;
}

std::string helpText(const char* executable_name)
{
    const std::string name = executable_name == nullptr ? "frehg2" : executable_name;
    return "Usage: " + name + " <config.yaml>\n"
           "\n"
           "Frehg2 Phase 1 infrastructure stub.\n"
           "Later phases will parse YAML, allocate domains, and run solvers.\n";
}

}  // namespace frehg2
