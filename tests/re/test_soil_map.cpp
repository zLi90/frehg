#include "re/ReSolver.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    int result = 0;
    frehg2::GridSpec spec;
    spec.nx = 2;
    spec.ny = 2;
    spec.nz = 3;
    spec.dx = 1.0;
    spec.dy = 1.0;
    spec.dz = 0.1;
    const frehg2::Grid grid(spec);

    const std::string filename = "phase10_soil_map.asc";
    {
        std::ofstream out(filename);
        out << "ncols 2\n";
        out << "nrows 2\n";
        out << "xllcorner 0.0\n";
        out << "yllcorner 0.0\n";
        out << "cellsize 1.0\n";
        out << "NODATA_value -9999\n";
        out << "2 3\n";
        out << "0 1\n";
    }

    const auto surface_soil = frehg2::ReSolver::readSurfaceSoilMap(grid, filename);
    if (surface_soil[grid.getSurfaceIndex(0, 0)] != 0 ||
        surface_soil[grid.getSurfaceIndex(1, 0)] != 1 ||
        surface_soil[grid.getSurfaceIndex(0, 1)] != 2 ||
        surface_soil[grid.getSurfaceIndex(1, 1)] != 3) {
        std::cerr << "surface soil map orientation mismatch\n";
        result = 1;
    }

    const auto soil_id = frehg2::ReSolver::expandSurfaceSoilMap(grid, surface_soil);
    if (soil_id[grid.getIndex(1, 1, 0)] != 3 ||
        soil_id[grid.getIndex(1, 1, 2)] != 3 ||
        soil_id[grid.getIndex(0, 0, 1)] != 0) {
        std::cerr << "surface soil IDs were not expanded through groundwater columns\n";
        result = 1;
    }

    const std::string volume_filename = "phase10_soil_map_3d.asc";
    {
        std::ofstream out(volume_filename);
        out << "ncols 2\n";
        out << "nrows 2\n";
        out << "nlayers 3\n";
        out << "xllcorner 0.0\n";
        out << "yllcorner 0.0\n";
        out << "cellsize 1.0\n";
        out << "NODATA_value -9999\n";
        out << "2 3\n";
        out << "0 1\n";
        out << "6 7\n";
        out << "4 5\n";
        out << "10 11\n";
        out << "8 9\n";
    }

    const auto volume_soil = frehg2::ReSolver::readSoilMap(grid, volume_filename);
    if (volume_soil[grid.getIndex(0, 0, 0)] != 0 ||
        volume_soil[grid.getIndex(0, 0, 1)] != 4 ||
        volume_soil[grid.getIndex(0, 0, 2)] != 8 ||
        volume_soil[grid.getIndex(1, 1, 1)] != 7) {
        std::cerr << "3D soil map orientation mismatch\n";
        result = 1;
    }

    const std::string b3_soil_map = std::string(FREHG2_SOURCE_DIR) +
                                    "/benchmarks/b3-kirkland/rasters/soilID.asc";
    frehg2::GridSpec b3_spec;
    b3_spec.nx = 50;
    b3_spec.ny = 1;
    b3_spec.nz = 30;
    b3_spec.dx = 0.1;
    b3_spec.dy = 0.1;
    b3_spec.dz = 0.1;
    const frehg2::Grid b3_grid(b3_spec);
    const auto b3_soil = frehg2::ReSolver::readSoilMap(b3_grid, b3_soil_map);
    if (b3_soil[b3_grid.getIndex(25, 0, 0)] != 0 ||
        b3_soil[b3_grid.getIndex(25, 0, 10)] != 1 ||
        b3_soil[b3_grid.getIndex(25, 0, 20)] != 0) {
        std::cerr << "b3 3D layered soil map was not loaded correctly\n";
        result = 1;
    }

    frehg2::GridSpec lateral_spec;
    lateral_spec.nx = 2;
    lateral_spec.ny = 1;
    lateral_spec.nz = 1;
    lateral_spec.dx = 1.0;
    lateral_spec.dy = 1.0;
    lateral_spec.dz = 1.0;
    frehg2::ReParameters lateral_params;
    lateral_params.ksx = 1.0;
    lateral_params.ksy = 1.0;
    lateral_params.ksz = 1.0;
    lateral_params.soil_table = {
        {10.0, 10.0, 10.0, 1.0, 2.0, 0.4, 0.0, 0.0},
        {2.0, 2.0, 2.0, 1.0, 2.0, 0.4, 0.0, 0.0},
    };
    lateral_params.soil_id = {0, 1};
    lateral_params.init_h = 1.0;
    lateral_params.init_wc = -1.0;
    frehg2::ReSolver lateral_solver(frehg2::Grid(lateral_spec), lateral_params);
    auto lateral_state = lateral_solver.initializeLegacyState(0.0, -1.0);
    lateral_solver.computeConductivityFaces(lateral_state);
    if (lateral_state.kx[lateral_solver.grid().getIndex(0, 0, 0)] != 6.0) {
        std::cerr << "lateral face conductivity did not use adjacent soil types\n";
        result = 1;
    }

    frehg2::GridSpec column_spec;
    column_spec.nx = 1;
    column_spec.ny = 1;
    column_spec.nz = 2;
    column_spec.dx = 1.0;
    column_spec.dy = 1.0;
    column_spec.dz = 1.0;
    frehg2::ReParameters params;
    params.ksx = 1.0;
    params.ksy = 1.0;
    params.ksz = 1.0;
    params.soil_table = {
        {1.0, 1.0, 10.0, 1.0, 2.0, 0.4, 0.0, 0.0},
        {1.0, 1.0, 2.0, 1.0, 2.0, 0.4, 0.0, 0.0},
    };
    params.soil_id = {0, 1};
    params.init_h = 1.0;
    params.init_wc = -1.0;
    frehg2::ReSolver flux_solver(frehg2::Grid(column_spec), params);
    auto flux_state = flux_solver.initializeLegacyState(0.0, -2.0);
    flux_state.h[0] = 1.0;
    flux_state.h[1] = 1.0;
    flux_solver.computeConductivityFaces(flux_state);
    flux_solver.computeDarcyFlux(flux_state);
    if (flux_state.kz[0] != 6.0 || flux_state.qz[0] != -6.0) {
        std::cerr << "vertical Darcy flux did not use heterogeneous face conductivity\n";
        result = 1;
    }

#ifdef USE_KOKKOS
    Kokkos::finalize();
#endif
    return result;
}
