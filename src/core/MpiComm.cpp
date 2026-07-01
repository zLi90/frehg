#include "core/MpiComm.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace frehg2 {

MpiComm::MpiComm(int global_nx, int global_ny, int mpi_nx, int mpi_ny, MPI_Comm comm)
    : comm_(comm),
      global_nx_(global_nx),
      global_ny_(global_ny),
      mpi_nx_(mpi_nx),
      mpi_ny_(mpi_ny) {
  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &size_);
  if (mpi_nx < 1 || mpi_ny < 1) {
    throw std::runtime_error("MpiComm: mpi_nx and mpi_ny must be >= 1");
  }
  if (mpi_nx * mpi_ny != size_) {
    throw std::runtime_error(
        "MpiComm: mpi_nx*mpi_ny (" + std::to_string(mpi_nx * mpi_ny) +
        ") must equal the communicator size (" + std::to_string(size_) + ")");
  }
  if (mpi_nx > global_nx || mpi_ny > global_ny) {
    throw std::runtime_error("MpiComm: process grid finer than the cell grid");
  }
}

int MpiComm::partSize(int part, int n, int P) {
  const int base = n / P;
  const int rem = n % P;
  return base + (part < rem ? 1 : 0);
}

int MpiComm::partStart(int part, int n, int P) {
  const int base = n / P;
  const int rem = n % P;
  return part * base + std::min(part, rem);
}

int MpiComm::partOf(int global_index, int n, int P) {
  const int base = n / P;
  const int rem = n % P;
  const int cut = rem * (base + 1);
  if (global_index < cut) {
    return global_index / (base + 1);
  }
  return rem + (global_index - cut) / base;
}

std::tuple<int, int, int> MpiComm::globalToLocal(int gi, int gj) const {
  const int px = partOf(gi, global_nx_, mpi_nx_);
  const int py = partOf(gj, global_ny_, mpi_ny_);
  const int owner = py * mpi_nx_ + px;
  const int li = gi - partStart(px, global_nx_, mpi_nx_);
  const int lj = gj - partStart(py, global_ny_, mpi_ny_);
  return {li, lj, owner};
}

std::tuple<int, int> MpiComm::localToGlobal(int li, int lj) const {
  return {i0() + li, j0() + lj};
}

void MpiComm::haloExchange2D(const RealArr1D& field, int nx_with_halo, int ny_with_halo,
                             HaloTags tags) const {
  const int lnx = localNx();
  const int lny = localNy();
  if (nx_with_halo != lnx + 2 || ny_with_halo != lny + 2) {
    throw std::runtime_error("haloExchange2D: field halo extents do not match this rank");
  }

  auto h = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(h, field);
  real* f = h.data();
  auto at = [&](int c, int r) -> real& { return f[c + r * nx_with_halo]; };

  // ---- x-direction (left/right) ----
  std::vector<real> send_l(lny), send_r(lny), recv_l(lny), recv_r(lny);
  for (int r = 0; r < lny; ++r) {
    send_l[r] = at(1, r + 1);        // leftmost interior column
    send_r[r] = at(lnx, r + 1);      // rightmost interior column
  }
  MPI_Sendrecv(send_l.data(), lny, MPI_DOUBLE, leftRank(), tags.left, recv_r.data(), lny,
               MPI_DOUBLE, rightRank(), tags.left, comm_, MPI_STATUS_IGNORE);
  MPI_Sendrecv(send_r.data(), lny, MPI_DOUBLE, rightRank(), tags.right, recv_l.data(), lny,
               MPI_DOUBLE, leftRank(), tags.right, comm_, MPI_STATUS_IGNORE);
  if (leftRank() != MPI_PROC_NULL) {
    for (int r = 0; r < lny; ++r) at(0, r + 1) = recv_l[r];
  }
  if (rightRank() != MPI_PROC_NULL) {
    for (int r = 0; r < lny; ++r) at(lnx + 1, r + 1) = recv_r[r];
  }

  // ---- y-direction (down/up): include ghost columns so corners are filled ----
  const int width = nx_with_halo;
  std::vector<real> send_d(width), send_u(width), recv_d(width), recv_u(width);
  for (int c = 0; c < width; ++c) {
    send_d[c] = at(c, 1);            // bottom interior row (with x-ghosts already filled)
    send_u[c] = at(c, lny);          // top interior row
  }
  MPI_Sendrecv(send_d.data(), width, MPI_DOUBLE, downRank(), tags.down, recv_u.data(),
               width, MPI_DOUBLE, upRank(), tags.down, comm_, MPI_STATUS_IGNORE);
  MPI_Sendrecv(send_u.data(), width, MPI_DOUBLE, upRank(), tags.up, recv_d.data(), width,
               MPI_DOUBLE, downRank(), tags.up, comm_, MPI_STATUS_IGNORE);
  if (downRank() != MPI_PROC_NULL) {
    for (int c = 0; c < width; ++c) at(c, 0) = recv_d[c];
  }
  if (upRank() != MPI_PROC_NULL) {
    for (int c = 0; c < width; ++c) at(c, lny + 1) = recv_u[c];
  }

  Kokkos::deep_copy(field, h);
}

void MpiComm::haloExchange3D(const RealArr1D& field, int nx_with_halo, int ny_with_halo,
                             int nz, HaloTags tags) const {
  const int plane = nx_with_halo * ny_with_halo;
  for (int k = 0; k < nz; ++k) {
    // Interior z layer k lives at storage offset (k+1)*plane (k=-1/nz are BC ghosts).
    RealArr1D layer(field.data() + static_cast<size_t>(k + 1) * plane, plane);
    haloExchange2D(layer, nx_with_halo, ny_with_halo, tags);
  }
}

void MpiComm::gatherToRank0(const RealArr1D& local, int nx_with_halo, int ny_with_halo,
                            const RealArr1DHost& global) const {
  const int lnx = localNx();
  const int lny = localNy();
  if (nx_with_halo != lnx + 2 || ny_with_halo != lny + 2) {
    throw std::runtime_error("gatherToRank0: field halo extents do not match this rank");
  }

  auto h = Kokkos::create_mirror_view(local);
  Kokkos::deep_copy(h, local);
  const real* f = h.data();

  // Pack this rank's interior block (row-major within block).
  std::vector<real> block(static_cast<size_t>(lnx) * lny);
  for (int r = 0; r < lny; ++r) {
    for (int c = 0; c < lnx; ++c) {
      block[c + r * lnx] = f[(c + 1) + (r + 1) * nx_with_halo];
    }
  }
  const int header[4] = {i0(), j0(), lnx, lny};

  if (rank_ != 0) {
    MPI_Send(header, 4, MPI_INT, 0, 9001, comm_);
    MPI_Send(block.data(), static_cast<int>(block.size()), MPI_DOUBLE, 0, 9002, comm_);
    return;
  }

  // Rank 0 places every rank's block into the global array.
  auto place = [&](const int hdr[4], const std::vector<real>& b) {
    const int bi0 = hdr[0], bj0 = hdr[1], bnx = hdr[2], bny = hdr[3];
    for (int r = 0; r < bny; ++r) {
      for (int c = 0; c < bnx; ++c) {
        const int gi = bi0 + c;
        const int gj = bj0 + r;
        global(gi + gj * global_nx_) = b[c + r * bnx];
      }
    }
  };
  place(header, block);
  for (int src = 1; src < size_; ++src) {
    int hdr[4];
    MPI_Recv(hdr, 4, MPI_INT, src, 9001, comm_, MPI_STATUS_IGNORE);
    std::vector<real> b(static_cast<size_t>(hdr[2]) * hdr[3]);
    MPI_Recv(b.data(), static_cast<int>(b.size()), MPI_DOUBLE, src, 9002, comm_,
             MPI_STATUS_IGNORE);
    place(hdr, b);
  }
}

}  // namespace frehg2
