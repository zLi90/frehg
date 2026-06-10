#include "core/MpiComm.hpp"

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdexcept>
#include <vector>

namespace frehg2 {

namespace {

int mpiProcNull()
{
#ifdef USE_MPI
    return MPI_PROC_NULL;
#else
    return -1;
#endif
}

}  // namespace

MpiComm::MpiComm(int mpi_nx, int mpi_ny)
    : mpi_nx_(mpi_nx),
      mpi_ny_(mpi_ny)
{
    if (mpi_nx_ <= 0 || mpi_ny_ <= 0) {
        throw std::invalid_argument("MPI decomposition dimensions must be positive");
    }

#ifdef USE_MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized != 0) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
    }
#endif

    if (size_ != mpi_nx_ * mpi_ny_) {
        throw std::invalid_argument("MPI decomposition does not match communicator size");
    }
}

int MpiComm::rank() const noexcept
{
    return rank_;
}

int MpiComm::size() const noexcept
{
    return size_;
}

int MpiComm::mpiNx() const noexcept
{
    return mpi_nx_;
}

int MpiComm::mpiNy() const noexcept
{
    return mpi_ny_;
}

int MpiComm::rankX() const noexcept
{
    return rank_ % mpi_nx_;
}

int MpiComm::rankY() const noexcept
{
    return rank_ / mpi_nx_;
}

int MpiComm::neighbor(Direction direction) const noexcept
{
    switch (direction) {
    case Direction::XM:
        return rankX() == 0 ? mpiProcNull() : rank_ - 1;
    case Direction::XP:
        return rankX() == mpi_nx_ - 1 ? mpiProcNull() : rank_ + 1;
    case Direction::YM:
        return rankY() == 0 ? mpiProcNull() : rank_ - mpi_nx_;
    case Direction::YP:
        return rankY() == mpi_ny_ - 1 ? mpiProcNull() : rank_ + mpi_nx_;
    case Direction::ZP:
    case Direction::ZM:
        return mpiProcNull();
    }

    return mpiProcNull();
}

#ifdef USE_KOKKOS
void MpiComm::exchangeHalos2D(realDeviceArr values, const Grid& grid) const
{
    if (values.extent(0) < grid.nSurfaceCellMem()) {
        throw std::invalid_argument("2D halo exchange view is too small");
    }
    if (size_ == 1) {
        return;
    }

    auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);

#ifdef USE_MPI
    std::vector<real> send_xp(static_cast<std::size_t>(grid.ny()));
    std::vector<real> send_xm(static_cast<std::size_t>(grid.ny()));
    std::vector<real> recv_xp(static_cast<std::size_t>(grid.ny()));
    std::vector<real> recv_xm(static_cast<std::size_t>(grid.ny()));
    for (int j = 0; j < grid.ny(); ++j) {
        send_xp[static_cast<std::size_t>(j)] = host(grid.getSurfaceIndex(grid.nx() - 1, j));
        send_xm[static_cast<std::size_t>(j)] = host(grid.getSurfaceIndex(0, j));
    }

    MPI_Sendrecv(send_xp.data(), grid.ny(), MPI_DOUBLE, neighbor(Direction::XP), 101,
                 recv_xm.data(), grid.ny(), MPI_DOUBLE, neighbor(Direction::XM), 101,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_xm.data(), grid.ny(), MPI_DOUBLE, neighbor(Direction::XM), 102,
                 recv_xp.data(), grid.ny(), MPI_DOUBLE, neighbor(Direction::XP), 102,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int j = 0; j < grid.ny(); ++j) {
        if (neighbor(Direction::XP) != MPI_PROC_NULL) {
            host(grid.surfaceGhostIndex(Direction::XP, grid.nx() - 1, j)) =
                recv_xp[static_cast<std::size_t>(j)];
        }
        if (neighbor(Direction::XM) != MPI_PROC_NULL) {
            host(grid.surfaceGhostIndex(Direction::XM, 0, j)) =
                recv_xm[static_cast<std::size_t>(j)];
        }
    }

    std::vector<real> send_yp(static_cast<std::size_t>(grid.nx()));
    std::vector<real> send_ym(static_cast<std::size_t>(grid.nx()));
    std::vector<real> recv_yp(static_cast<std::size_t>(grid.nx()));
    std::vector<real> recv_ym(static_cast<std::size_t>(grid.nx()));
    for (int i = 0; i < grid.nx(); ++i) {
        send_yp[static_cast<std::size_t>(i)] = host(grid.getSurfaceIndex(i, grid.ny() - 1));
        send_ym[static_cast<std::size_t>(i)] = host(grid.getSurfaceIndex(i, 0));
    }

    MPI_Sendrecv(send_yp.data(), grid.nx(), MPI_DOUBLE, neighbor(Direction::YP), 201,
                 recv_ym.data(), grid.nx(), MPI_DOUBLE, neighbor(Direction::YM), 201,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_ym.data(), grid.nx(), MPI_DOUBLE, neighbor(Direction::YM), 202,
                 recv_yp.data(), grid.nx(), MPI_DOUBLE, neighbor(Direction::YP), 202,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < grid.nx(); ++i) {
        if (neighbor(Direction::YP) != MPI_PROC_NULL) {
            host(grid.surfaceGhostIndex(Direction::YP, i, grid.ny() - 1)) =
                recv_yp[static_cast<std::size_t>(i)];
        }
        if (neighbor(Direction::YM) != MPI_PROC_NULL) {
            host(grid.surfaceGhostIndex(Direction::YM, i, 0)) =
                recv_ym[static_cast<std::size_t>(i)];
        }
    }
#endif

    Kokkos::deep_copy(values, host);
}

void MpiComm::exchangeHalos3D(realDeviceArr values, const Grid& grid) const
{
    if (values.extent(0) < grid.nCellMem()) {
        throw std::invalid_argument("3D halo exchange view is too small");
    }
    if (size_ == 1) {
        return;
    }

    auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), values);

#ifdef USE_MPI
    const int yz_count = grid.ny() * grid.nz();
    std::vector<real> send_xp(static_cast<std::size_t>(yz_count));
    std::vector<real> send_xm(static_cast<std::size_t>(yz_count));
    std::vector<real> recv_xp(static_cast<std::size_t>(yz_count));
    std::vector<real> recv_xm(static_cast<std::size_t>(yz_count));
    for (int j = 0; j < grid.ny(); ++j) {
        for (int k = 0; k < grid.nz(); ++k) {
            const auto offset = static_cast<std::size_t>(j * grid.nz() + k);
            send_xp[offset] = host(grid.getIndex(grid.nx() - 1, j, k));
            send_xm[offset] = host(grid.getIndex(0, j, k));
        }
    }

    MPI_Sendrecv(send_xp.data(), yz_count, MPI_DOUBLE, neighbor(Direction::XP), 301,
                 recv_xm.data(), yz_count, MPI_DOUBLE, neighbor(Direction::XM), 301,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_xm.data(), yz_count, MPI_DOUBLE, neighbor(Direction::XM), 302,
                 recv_xp.data(), yz_count, MPI_DOUBLE, neighbor(Direction::XP), 302,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int j = 0; j < grid.ny(); ++j) {
        for (int k = 0; k < grid.nz(); ++k) {
            const auto offset = static_cast<std::size_t>(j * grid.nz() + k);
            if (neighbor(Direction::XP) != MPI_PROC_NULL) {
                host(grid.groundwaterGhostIndex(Direction::XP, grid.nx() - 1, j, k)) =
                    recv_xp[offset];
            }
            if (neighbor(Direction::XM) != MPI_PROC_NULL) {
                host(grid.groundwaterGhostIndex(Direction::XM, 0, j, k)) = recv_xm[offset];
            }
        }
    }

    const int xz_count = grid.nx() * grid.nz();
    std::vector<real> send_yp(static_cast<std::size_t>(xz_count));
    std::vector<real> send_ym(static_cast<std::size_t>(xz_count));
    std::vector<real> recv_yp(static_cast<std::size_t>(xz_count));
    std::vector<real> recv_ym(static_cast<std::size_t>(xz_count));
    for (int i = 0; i < grid.nx(); ++i) {
        for (int k = 0; k < grid.nz(); ++k) {
            const auto offset = static_cast<std::size_t>(i * grid.nz() + k);
            send_yp[offset] = host(grid.getIndex(i, grid.ny() - 1, k));
            send_ym[offset] = host(grid.getIndex(i, 0, k));
        }
    }

    MPI_Sendrecv(send_yp.data(), xz_count, MPI_DOUBLE, neighbor(Direction::YP), 401,
                 recv_ym.data(), xz_count, MPI_DOUBLE, neighbor(Direction::YM), 401,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_ym.data(), xz_count, MPI_DOUBLE, neighbor(Direction::YM), 402,
                 recv_yp.data(), xz_count, MPI_DOUBLE, neighbor(Direction::YP), 402,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < grid.nx(); ++i) {
        for (int k = 0; k < grid.nz(); ++k) {
            const auto offset = static_cast<std::size_t>(i * grid.nz() + k);
            if (neighbor(Direction::YP) != MPI_PROC_NULL) {
                host(grid.groundwaterGhostIndex(Direction::YP, i, grid.ny() - 1, k)) =
                    recv_yp[offset];
            }
            if (neighbor(Direction::YM) != MPI_PROC_NULL) {
                host(grid.groundwaterGhostIndex(Direction::YM, i, 0, k)) = recv_ym[offset];
            }
        }
    }
#endif

    Kokkos::deep_copy(values, host);
}
#endif

}  // namespace frehg2
