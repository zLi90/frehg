#ifndef FREHG2_CORE_TYPES_HPP
#define FREHG2_CORE_TYPES_HPP

#include "core/define.hpp"

#include <array>
#include <string>

namespace frehg2 {

struct GridSpec {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    real dx = 1.0;
    real dy = 1.0;
    real dz = 1.0;
    real dz_multiplier = 1.0;
};

struct TimeSpec {
    real start = 0.0;
    real end = 0.0;
    real dt = 0.0;
    real output_interval = 0.0;
};

struct BoundaryCodes {
    std::array<int, 4> surface_water = {0, 0, 0, 0};
    std::array<int, 6> groundwater = {0, 0, 0, 0, 0, 0};
};

struct SimulationConfig {
    std::string input_path;
    SimMode mode = SimMode::COUPLED;
    GridSpec grid;
    TimeSpec time;
    BoundaryCodes boundary_codes;
};

struct BuildInfo {
    bool kokkos_enabled = false;
    bool mpi_enabled = false;
    bool petsc_enabled = false;
    bool hdf5_enabled = false;
    bool yaml_cpp_enabled = false;
};

BuildInfo getBuildInfo();
const char* simModeName(SimMode mode);

}  // namespace frehg2

#endif  // FREHG2_CORE_TYPES_HPP
