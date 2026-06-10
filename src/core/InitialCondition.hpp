#ifndef FREHG2_CORE_INITIAL_CONDITION_HPP
#define FREHG2_CORE_INITIAL_CONDITION_HPP

#include "core/Grid.hpp"

#include <string>
#include <vector>

namespace frehg2 {

enum class InitialConditionSource {
    CONSTANT,
    ASCII_RASTER,
    HDF5_VECTOR
};

struct InitialConditionSpec {
    std::string variable;
    InitialConditionSource source = InitialConditionSource::CONSTANT;
    real value = 0.0;
    std::string filename;
    std::string dataset;
};

class InitialCondition {
public:
    static std::vector<real> loadSurfaceField(
        const Grid& grid,
        const InitialConditionSpec& spec);

    static std::vector<real> loadGroundwaterField(
        const Grid& grid,
        const InitialConditionSpec& spec);

    static std::vector<real> readAsciiGroundwaterField(
        const Grid& grid,
        const std::string& filename);

private:
    static std::vector<real> readAsciiSurface(const Grid& grid, const std::string& filename);
    static std::vector<real> readHdf5Vector(const std::string& filename, const std::string& dataset);
};

}  // namespace frehg2

#endif  // FREHG2_CORE_INITIAL_CONDITION_HPP
