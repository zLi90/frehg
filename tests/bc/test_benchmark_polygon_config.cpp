#include "bc/BoundaryCondition.hpp"
#include "io/Config.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace {

std::string pathJoin(const std::string& directory, const std::string& filename)
{
    if (filename.empty() || filename.front() == '/') {
        return filename;
    }
    return directory + "/" + filename;
}

frehg2::BoundaryConditionType bcType(const std::string& value)
{
    if (value == "zero_gradient") {
        return frehg2::BoundaryConditionType::ZERO_GRADIENT;
    }
    if (value == "groundwater_head") {
        return frehg2::BoundaryConditionType::GROUNDWATER_HEAD;
    }
    if (value == "fixed_water_level") {
        return frehg2::BoundaryConditionType::FIXED_WATER_LEVEL;
    }
    return frehg2::BoundaryConditionType::ZERO_GRADIENT;
}

frehg2::SourceSinkType sourceType(const std::string& value)
{
    if (value == "rainfall") {
        return frehg2::SourceSinkType::RAINFALL;
    }
    if (value == "evapotranspiration") {
        return frehg2::SourceSinkType::EVAPOTRANSPIRATION;
    }
    if (value == "inflow") {
        return frehg2::SourceSinkType::INFLOW;
    }
    return frehg2::SourceSinkType::PUMPING;
}

frehg2::PolygonBoundaryCondition loadBc(const frehg2::Config& config, const std::string& path)
{
    return frehg2::PolygonBoundaryCondition{
        frehg2::Polygon::readFromFile(
            pathJoin(config.baseDirectory(), config.get<std::string>(path + ".polygon_file"))),
        bcType(config.get<std::string>(path + ".type")),
        config.getOr<frehg2::real>(path + ".value", 0.0),
    };
}

frehg2::PolygonSourceSink loadSource(const frehg2::Config& config, const std::string& path)
{
    return frehg2::PolygonSourceSink{
        frehg2::Polygon::readFromFile(
            pathJoin(config.baseDirectory(), config.get<std::string>(path + ".polygon_file"))),
        sourceType(config.get<std::string>(path + ".type")),
        config.getOr<frehg2::real>(path + ".rate", 0.0),
    };
}

}  // namespace

int main()
{
    const std::string root = FREHG2_SOURCE_DIR;

    const frehg2::Config sw(root + "/benchmarks/b1-sw/b1-sw.yaml");
    if (!sw.getOr<bool>("polygon_boundary_conditions.enable", false) ||
        !sw.getOr<bool>("polygon_sources.enable", false)) {
        std::cerr << "b1 polygon sections are not activated\n";
        return 1;
    }
    if (sw.indexedCount("polygon_boundary_conditions.surface") != 1 ||
        sw.indexedCount("polygon_sources.surface") != 1) {
        std::cerr << "b1 polygon section counts are wrong\n";
        return 1;
    }

    auto sw_spec = sw.gridSpec();
    sw_spec.nz = 1;
    const frehg2::Grid sw_grid(sw_spec);
    const auto sw_bc = loadBc(sw, "polygon_boundary_conditions.surface.0");
    const auto markers = frehg2::BoundaryCondition::markCells(sw_grid, {sw_bc});
    if (markers.size() != static_cast<std::size_t>(sw_grid.nSurfaceCell()) || markers.front() == 0) {
        std::cerr << "b1 polygon BC did not select the benchmark domain\n";
        return 1;
    }

    std::vector<frehg2::real> depth(static_cast<std::size_t>(sw_grid.nSurfaceCell()), 0.5);
    frehg2::BoundaryCondition::applySourceSink(
        sw_grid,
        {loadSource(sw, "polygon_sources.surface.0")},
        depth,
        sw.timeSpec().dt);
    for (const auto value : depth) {
        if (value != 0.5) {
            std::cerr << "b1 zero-rate polygon rainfall changed legacy-compatible depth\n";
            return 1;
        }
    }

    const frehg2::Config gw(root + "/benchmarks/b2-gw/b2-gw.yaml");
    if (!gw.getOr<bool>("polygon_boundary_conditions.enable", false) ||
        gw.indexedCount("polygon_boundary_conditions.groundwater") != 1) {
        std::cerr << "b2 groundwater polygon BC is not activated\n";
        return 1;
    }

    auto gw_spec = gw.gridSpec();
    gw_spec.nz = 4;
    const frehg2::Grid gw_grid(gw_spec);
    std::vector<frehg2::real> head(static_cast<std::size_t>(gw_grid.nCell()), -1.0);
    frehg2::BoundaryCondition::applyGroundwaterHead(
        gw_grid,
        {loadBc(gw, "polygon_boundary_conditions.groundwater.0")},
        head);
    for (const auto value : head) {
        if (value != gw.get<frehg2::real>("groundwater.htop")) {
            std::cerr << "b2 groundwater-head polygon did not cover the benchmark domain\n";
            return 1;
        }
    }

    return 0;
}
