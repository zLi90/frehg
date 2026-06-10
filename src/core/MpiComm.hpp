#ifndef FREHG2_CORE_MPI_COMM_HPP
#define FREHG2_CORE_MPI_COMM_HPP

#include "core/Grid.hpp"

namespace frehg2 {

class MpiComm {
public:
    MpiComm(int mpi_nx = 1, int mpi_ny = 1);

    int rank() const noexcept;
    int size() const noexcept;
    int mpiNx() const noexcept;
    int mpiNy() const noexcept;
    int rankX() const noexcept;
    int rankY() const noexcept;

    int neighbor(Direction direction) const noexcept;

#ifdef USE_KOKKOS
    void exchangeHalos2D(realDeviceArr values, const Grid& grid) const;
    void exchangeHalos3D(realDeviceArr values, const Grid& grid) const;
#endif

private:
    int rank_ = 0;
    int size_ = 1;
    int mpi_nx_ = 1;
    int mpi_ny_ = 1;
};

}  // namespace frehg2

#endif  // FREHG2_CORE_MPI_COMM_HPP
