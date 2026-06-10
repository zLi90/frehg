#include "io/TimeSeries.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>

int main()
{
    const auto path = std::filesystem::temp_directory_path() / "frehg2_timeseries_test.txt";

    {
        std::ofstream output(path);
        output << "# t value\n";
        output << "0, 1\n";
        output << "10 3\n";
        output << "20, 7\n";
    }

    try {
        frehg2::TimeSeries series;
        series.read(path.string());

        if (series.size() != 3) {
            return 1;
        }
        if (series.getValueAt(-1.0) != 1.0 || series.getValueAt(21.0) != 7.0) {
            return 1;
        }
        if (std::abs(series.getValueAt(5.0) - 2.0) > 1.0e-12) {
            return 1;
        }
        if (std::abs(series.getValueAt(15.0) - 5.0) > 1.0e-12) {
            return 1;
        }
    } catch (...) {
        return 1;
    }

    std::filesystem::remove(path);
    return 0;
}
