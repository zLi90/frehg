#ifndef FREHG2_CORE_DEFINE_HPP
#define FREHG2_CORE_DEFINE_HPP

#include <cstddef>
#include <string>

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

namespace frehg2 {

using real = double;
using index_t = std::size_t;

#ifdef USE_KOKKOS
using realArr = Kokkos::View<real*, Kokkos::HostSpace>;
using realArr2 = Kokkos::View<real**, Kokkos::HostSpace>;
using intArr = Kokkos::View<int*, Kokkos::HostSpace>;
using realDeviceArr = Kokkos::View<real*>;
using realDeviceArr2 = Kokkos::View<real**>;
using intDeviceArr = Kokkos::View<int*>;
#endif

enum class SimMode {
    SW_ONLY,
    GW_ONLY,
    COUPLED,
    SOLUTE
};

enum class BCType {
    DIRICHLET,
    NEUMANN,
    FREEFLOW,
    ZEROGRADIENT
};

}  // namespace frehg2

#endif  // FREHG2_CORE_DEFINE_HPP
