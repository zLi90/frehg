#ifndef FREHG2_IO_TIME_SERIES_HPP
#define FREHG2_IO_TIME_SERIES_HPP

#include "core/define.hpp"

#include <string>
#include <vector>

namespace frehg2 {

class TimeSeries {
public:
    void read(const std::string& filename);
    void addPoint(real time, real value);

    index_t size() const noexcept;
    real time(index_t index) const;
    real value(index_t index) const;
    real getValueAt(real query_time) const;

private:
    std::vector<real> times_;
    std::vector<real> values_;

    void ensureSorted() const;
};

}  // namespace frehg2

#endif  // FREHG2_IO_TIME_SERIES_HPP
