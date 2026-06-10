#include "io/TimeSeries.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace frehg2 {

namespace {

std::string normalizeDelimiters(std::string line)
{
    for (auto& ch : line) {
        if (ch == ',') {
            ch = ' ';
        }
    }
    return line;
}

}  // namespace

void TimeSeries::read(const std::string& filename)
{
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open time series: " + filename);
    }

    times_.clear();
    values_.clear();

    std::string line;
    while (std::getline(input, line)) {
        const auto comment = line.find('#');
        if (comment != std::string::npos) {
            line = line.substr(0, comment);
        }
        line = normalizeDelimiters(line);

        std::istringstream row(line);
        real t = 0.0;
        real v = 0.0;
        if (row >> t >> v) {
            addPoint(t, v);
        }
    }

    if (times_.empty()) {
        throw std::runtime_error("time series contains no data");
    }
    ensureSorted();
}

void TimeSeries::addPoint(real time, real value)
{
    times_.push_back(time);
    values_.push_back(value);
}

index_t TimeSeries::size() const noexcept
{
    return times_.size();
}

real TimeSeries::time(index_t index) const
{
    if (index >= times_.size()) {
        throw std::out_of_range("time-series time index out of range");
    }
    return times_[index];
}

real TimeSeries::value(index_t index) const
{
    if (index >= values_.size()) {
        throw std::out_of_range("time-series value index out of range");
    }
    return values_[index];
}

real TimeSeries::getValueAt(real query_time) const
{
    if (times_.empty()) {
        throw std::runtime_error("cannot interpolate an empty time series");
    }
    ensureSorted();

    if (query_time <= times_.front()) {
        return values_.front();
    }
    if (query_time >= times_.back()) {
        return values_.back();
    }

    const auto upper = std::upper_bound(times_.begin(), times_.end(), query_time);
    const auto index = static_cast<std::size_t>(upper - times_.begin());
    const real t0 = times_[index - 1];
    const real t1 = times_[index];
    const real v0 = values_[index - 1];
    const real v1 = values_[index];

    const real fraction = (query_time - t0) / (t1 - t0);
    return v0 + fraction * (v1 - v0);
}

void TimeSeries::ensureSorted() const
{
    if (!std::is_sorted(times_.begin(), times_.end())) {
        throw std::runtime_error("time-series times must be sorted");
    }
}

}  // namespace frehg2
