#ifndef FREHG2_CORE_MONITOR_HPP
#define FREHG2_CORE_MONITOR_HPP

#include "core/Grid.hpp"

#include <string>
#include <vector>

namespace frehg2 {

struct MonitorStats {
    real sum = 0.0;
    real min = 0.0;
    real max = 0.0;
    real mean = 0.0;
    std::size_t count = 0;
};

struct MassBalance {
    real old_mass = 0.0;
    real new_mass = 0.0;
    real boundary_flux = 0.0;
    real source_sink = 0.0;
    real dt = 0.0;

    real residual() const noexcept;
};

class Monitor {
public:
    static real pointValue(const std::vector<real>& values, index_t cell);
    static MonitorStats summarizeCells(
        const std::vector<real>& values,
        const std::vector<index_t>& cells);
    static MassBalance massBalance(
        real old_mass,
        real new_mass,
        real boundary_flux,
        real source_sink,
        real dt);
    static void writePolygonHdf5(
        const std::string& filename,
        const std::string& group_name,
        const std::vector<index_t>& cells,
        const std::vector<real>& values);
};

}  // namespace frehg2

#endif  // FREHG2_CORE_MONITOR_HPP
