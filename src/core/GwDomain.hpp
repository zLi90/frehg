#ifndef FREHG2_CORE_GW_DOMAIN_HPP
#define FREHG2_CORE_GW_DOMAIN_HPP

#include "core/Grid.hpp"

namespace frehg2 {

class GwDomain {
public:
    explicit GwDomain(Grid grid);

    const Grid& grid() const noexcept;

    int nzGlob() const noexcept;
    real dzMultiplier() const noexcept;
    index_t size() const noexcept;

private:
    Grid grid_;

public:
#ifdef USE_KOKKOS
    intDeviceArr soilID;
    realDeviceArr dz;
#endif

    void setUniformSoil(int soil_id);
};

}  // namespace frehg2

#endif  // FREHG2_CORE_GW_DOMAIN_HPP
