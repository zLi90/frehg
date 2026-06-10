#include "io/AsciiRaster.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#include <filesystem>
#include <fstream>

int main(int argc, char** argv)
{
#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#else
    (void)argc;
    (void)argv;
#endif

    const auto path = std::filesystem::temp_directory_path() / "frehg2_ascii_raster_test.asc";
    int result = 0;

    {
        std::ofstream output(path);
        output << "ncols 3\n";
        output << "nrows 2\n";
        output << "xllcorner 10\n";
        output << "yllcorner 20\n";
        output << "cellsize 5\n";
        output << "NODATA_value -9999\n";
        output << "1 2 -9999\n";
        output << "4 5 6\n";
    }

    try {
        frehg2::AsciiRaster raster;
        raster.read(path.string());

        if (raster.ncols() != 3 || raster.nrows() != 2 || raster.size() != 6) {
            result = 1;
        }
        if (raster.xllcorner() != 10.0 || raster.yllcorner() != 20.0 ||
            raster.cellsize() != 5.0) {
            result = 1;
        }

        if (raster.value(0) != 4.0 || raster.value(1) != 5.0 || raster.value(2) != 6.0) {
            result = 1;
        }
        if (raster.value(5) != -9999.0 || raster.active(5) != 0) {
            result = 1;
        }

#ifdef USE_KOKKOS
        if (raster.dataView().extent(0) != 6 || raster.activeMaskView().extent(0) != 6) {
            result = 1;
        }
#endif
    } catch (...) {
        result = 1;
    }

    std::filesystem::remove(path);

#ifdef USE_KOKKOS
    Kokkos::finalize();
#endif

    return result;
}
