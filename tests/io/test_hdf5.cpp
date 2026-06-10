#include "io/Hdf5Reader.hpp"
#include "io/Hdf5Writer.hpp"

#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <vector>

int main()
{
    const auto path = std::filesystem::temp_directory_path() / "frehg2_phase3_test.h5";
    frehg2::Hdf5Writer writer(path.string());
    frehg2::Hdf5Reader reader(path.string());

    if (writer.filename() != path.string() || reader.filename() != path.string()) {
        return 1;
    }

#ifdef USE_HDF5
    if (!writer.isEnabled() || !reader.isEnabled()) {
        return 1;
    }
    writer.writeMetadata("Phase 3 Test", "0.1.0");
    writer.writeVector("/surface/water_depth", std::vector<frehg2::real>{1.0, 2.0, 4.0});

    const auto values = reader.readVector("/surface/water_depth");
    if (values.size() != 3) {
        return 1;
    }
    if (std::abs(values[0] - 1.0) > 1.0e-12 ||
        std::abs(values[1] - 2.0) > 1.0e-12 ||
        std::abs(values[2] - 4.0) > 1.0e-12) {
        return 1;
    }
#else
    if (writer.isEnabled() || reader.isEnabled()) {
        return 1;
    }

    try {
        writer.writeVector("/surface/water_depth", std::vector<frehg2::real>{1.0, 2.0});
        return 1;
    } catch (const std::runtime_error&) {
    }

    try {
        (void)reader.readVector("/surface/water_depth");
        return 1;
    } catch (const std::runtime_error&) {
    }
#endif

    std::filesystem::remove(path);
    return 0;
}
