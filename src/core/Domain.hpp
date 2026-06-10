#ifndef FREHG2_CORE_DOMAIN_HPP
#define FREHG2_CORE_DOMAIN_HPP

#include "core/Grid.hpp"

namespace frehg2 {

class Domain {
public:
    explicit Domain(Grid grid);

    const Grid& grid() const noexcept;

private:
    Grid grid_;

public:
#ifdef USE_KOKKOS
    realDeviceArr z;
    realDeviceArr area;
    intDeviceArr actMask;
#endif

    index_t size() const noexcept;
    void setUniformBed(real bed_elevation);
    void setAllActive(int active = 1);
};

}  // namespace frehg2

#endif  // FREHG2_CORE_DOMAIN_HPP
