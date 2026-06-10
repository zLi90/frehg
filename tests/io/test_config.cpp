#include "io/Config.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main()
{
    const std::string root = FREHG2_SOURCE_DIR;

    try {
        const frehg2::Config sw(root + "/benchmarks/b1-sw/b1-sw.yaml");
        sw.validate();
        if (sw.get<std::string>("simulation.id") != "MAXP1") {
            return 1;
        }
        if (sw.get<int>("domain.nx") != 1 || sw.get<int>("domain.ny") != 10) {
            return 1;
        }
        if (!sw.get<bool>("surface_water.enable")) {
            return 1;
        }
        if (sw.get<bool>("groundwater.enable")) {
            return 1;
        }
        if (!sw.getOr<bool>("groundwater.async", true)) {
            return 1;
        }
        if (sw.get<std::vector<int>>("surface_water.bc_type").size() != 4) {
            return 1;
        }
        if (!sw.getOr<bool>("polygon_boundary_conditions.enable", false) ||
            sw.indexedCount("polygon_boundary_conditions.surface") != 1 ||
            sw.indexedCount("polygon_sources.surface") != 1) {
            return 1;
        }
        if (sw.get<std::string>("polygon_boundary_conditions.surface.0.polygon_file") !=
            "polygons/domain.poly") {
            return 1;
        }

        const auto sw_grid = sw.gridSpec();
        if (sw_grid.nx != 1 || sw_grid.ny != 10 || sw_grid.dx != 80.0) {
            return 1;
        }
        const auto sw_time = sw.timeSpec();
        if (sw_time.dt != 5.0 || sw_time.end != 18000.0) {
            return 1;
        }

        const frehg2::Config gw(root + "/benchmarks/b2-gw/b2-gw.yaml");
        gw.validate();
        if (gw.get<std::string>("simulation.id") != "warri") {
            return 1;
        }
        if (gw.get<bool>("surface_water.enable")) {
            return 1;
        }
        if (!gw.get<bool>("groundwater.enable")) {
            return 1;
        }
        if (!gw.getOr<bool>("groundwater.async", true)) {
            return 1;
        }
        const auto gw_bc = gw.get<std::vector<int>>("groundwater.bc_type");
        if (gw_bc.size() != 6 || gw_bc[4] != 0 || gw_bc[5] != 1) {
            return 1;
        }
        if (gw.get<frehg2::real>("groundwater.Ksz") <= 0.0) {
            return 1;
        }
        if (gw.get<std::string>("groundwater.soil_map_file") != "soil/soil_map_uniform.asc" ||
            gw.indexedCount("groundwater.soil_types") != 1 ||
            gw.get<int>("groundwater.soil_types.0.id") != 0 ||
            gw.get<frehg2::real>("groundwater.soil_types.0.Ksz") !=
                gw.get<frehg2::real>("groundwater.Ksz")) {
            return 1;
        }
        if (!gw.getOr<bool>("polygon_boundary_conditions.enable", false) ||
            gw.indexedCount("polygon_boundary_conditions.groundwater") != 1) {
            return 1;
        }
        if (gw.get<std::string>("polygon_boundary_conditions.groundwater.0.type") !=
            "groundwater_head") {
            return 1;
        }

        const frehg2::Config b3(root + "/benchmarks/b3-kirkland/b3-kirkland.yaml");
        b3.validate();
        if (b3.get<std::string>("initial_conditions.groundwater.hydraulic_head.source") !=
            "ascii_raster_3d") {
            return 1;
        }
        if (b3.indexedCount("boundary_conditions.groundwater") != 1 ||
            b3.get<std::string>("boundary_conditions.groundwater.0.type") != "fixed_flow_rate") {
            return 1;
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
