#ifndef FREHG2_CORE_GW_STATE_HPP
#define FREHG2_CORE_GW_STATE_HPP

#include "core/Grid.hpp"

namespace frehg2 {

class GwState {
public:
    explicit GwState(const Grid& grid);

    index_t size() const noexcept;
    void fill(real hydraulic_head, real water_content);
    void swapTimeLevels();

private:
    index_t size_ = 0;

public:
#ifdef USE_KOKKOS
    realDeviceArr h_old;
    realDeviceArr h_new;
    realDeviceArr wc_old;
    realDeviceArr wc_new;
    realDeviceArr wc_excess;
    realDeviceArr2 k;
    realDeviceArr2 q;
#endif
};

}  // namespace frehg2

#endif  // FREHG2_CORE_GW_STATE_HPP
