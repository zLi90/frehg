#ifndef FREHG2_CORE_STATE_HPP
#define FREHG2_CORE_STATE_HPP

#include "core/Grid.hpp"

namespace frehg2 {

class State {
public:
    explicit State(const Grid& grid);

    index_t size() const noexcept;
    void fill(real water_depth, real momentum_x, real momentum_y);
    void swapTimeLevels();

private:
    index_t size_ = 0;

public:
#ifdef USE_KOKKOS
    realDeviceArr h_old;
    realDeviceArr h_new;
    realDeviceArr hu_old;
    realDeviceArr hu_new;
    realDeviceArr hv_old;
    realDeviceArr hv_new;
    realDeviceArr z;
    realDeviceArr roughness;
    realDeviceArr qss;
#endif
};

}  // namespace frehg2

#endif  // FREHG2_CORE_STATE_HPP
