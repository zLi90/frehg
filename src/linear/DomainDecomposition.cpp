#include "linear/DomainDecomposition.hpp"

namespace frehg2 {

// ===================================== Decomp2D =====================================

Decomp2D::Decomp2D(const MpiComm& comm)
    : DecompBase(comm), nx_(comm.globalNx()), ny_(comm.globalNy()) {}

std::tuple<int, int, int, int> Decomp2D::localCorners() const {
  return {comm_.i0(), comm_.j0(), comm_.localNx(), comm_.localNy()};
}

int Decomp2D::rowStartOfRank(int px, int py) const {
  const int mpi_nx = comm_.mpiNx();
  const int mpi_ny = comm_.mpiNy();
  const int lny = MpiComm::partSize(py, ny_, mpi_ny);
  // Cells in all full process-rows below py, plus partial cells to the left in row py.
  return nx_ * MpiComm::partStart(py, ny_, mpi_ny) +
         lny * MpiComm::partStart(px, nx_, mpi_nx);
}

int Decomp2D::globalRow(int gi, int gj) const {
  const int mpi_nx = comm_.mpiNx();
  const int mpi_ny = comm_.mpiNy();
  const int px = MpiComm::partOf(gi, nx_, mpi_nx);
  const int py = MpiComm::partOf(gj, ny_, mpi_ny);
  const int i0 = MpiComm::partStart(px, nx_, mpi_nx);
  const int j0 = MpiComm::partStart(py, ny_, mpi_ny);
  const int lnx = MpiComm::partSize(px, nx_, mpi_nx);
  const int local = (gj - j0) * lnx + (gi - i0);
  return rowStartOfRank(px, py) + local;
}

void Decomp2D::stencilColumns(int gi, int gj, int* cols, int& ncols) const {
  ncols = 0;
  cols[ncols++] = globalRow(gi, gj);                       // center
  if (gi - 1 >= 0) cols[ncols++] = globalRow(gi - 1, gj);  // W
  if (gi + 1 < nx_) cols[ncols++] = globalRow(gi + 1, gj); // E
  if (gj - 1 >= 0) cols[ncols++] = globalRow(gi, gj - 1);  // S
  if (gj + 1 < ny_) cols[ncols++] = globalRow(gi, gj + 1); // N
}

std::pair<int, int> Decomp2D::ownershipRange() const {
  const int start = rowStartOfRank(comm_.px(), comm_.py());
  return {start, start + comm_.localNx() * comm_.localNy()};
}

void Decomp2D::ownedRowStencil(int local_owned, int& global_row, int* cols,
                               int& ncols) const {
  const int lnx = comm_.localNx();
  const int li = local_owned % lnx;
  const int lj = local_owned / lnx;
  const int gi = comm_.i0() + li;
  const int gj = comm_.j0() + lj;
  global_row = globalRow(gi, gj);
  stencilColumns(gi, gj, cols, ncols);
}

void Decomp2D::haloExchange(const RealArr1D& field) const {
  comm_.haloExchange2D(field, comm_.localNx() + 2, comm_.localNy() + 2);
}

void Decomp2D::packOwned(const RealArr1D& field_halo, const RealArr1D& owned) const {
  const int lnx = comm_.localNx();
  const int lny = comm_.localNy();
  const int nxh = lnx + 2;
  auto fh = Kokkos::create_mirror_view(field_halo);
  Kokkos::deep_copy(fh, field_halo);
  auto oh = Kokkos::create_mirror_view(owned);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      oh(lj * lnx + li) = fh((li + 1) + (lj + 1) * nxh);
  Kokkos::deep_copy(owned, oh);
}

void Decomp2D::unpackOwned(const RealArr1D& owned, const RealArr1D& field_halo) const {
  const int lnx = comm_.localNx();
  const int lny = comm_.localNy();
  const int nxh = lnx + 2;
  auto oh = Kokkos::create_mirror_view(owned);
  Kokkos::deep_copy(oh, owned);
  auto fh = Kokkos::create_mirror_view(field_halo);
  Kokkos::deep_copy(fh, field_halo);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      fh((li + 1) + (lj + 1) * nxh) = oh(lj * lnx + li);
  Kokkos::deep_copy(field_halo, fh);
}

// ===================================== Decomp3D =====================================

Decomp3D::Decomp3D(const MpiComm& comm, int nz)
    : DecompBase(comm), nx_(comm.globalNx()), ny_(comm.globalNy()), nz_(nz) {}

std::tuple<int, int, int, int, int, int> Decomp3D::localCorners() const {
  return {comm_.i0(), comm_.j0(), 0, comm_.localNx(), comm_.localNy(), nz_};
}

int Decomp3D::rowStartOfRank(int px, int py) const {
  const int mpi_nx = comm_.mpiNx();
  const int mpi_ny = comm_.mpiNy();
  const int lny = MpiComm::partSize(py, ny_, mpi_ny);
  const int cells_2d = nx_ * MpiComm::partStart(py, ny_, mpi_ny) +
                       lny * MpiComm::partStart(px, nx_, mpi_nx);
  return cells_2d * nz_;
}

int Decomp3D::globalRow(int gi, int gj, int k) const {
  const int mpi_nx = comm_.mpiNx();
  const int mpi_ny = comm_.mpiNy();
  const int px = MpiComm::partOf(gi, nx_, mpi_nx);
  const int py = MpiComm::partOf(gj, ny_, mpi_ny);
  const int i0 = MpiComm::partStart(px, nx_, mpi_nx);
  const int j0 = MpiComm::partStart(py, ny_, mpi_ny);
  const int lnx = MpiComm::partSize(px, nx_, mpi_nx);
  const int local = ((gj - j0) * lnx + (gi - i0)) * nz_ + k;
  return rowStartOfRank(px, py) + local;
}

void Decomp3D::stencilColumns(int gi, int gj, int k, int* cols, int& ncols) const {
  ncols = 0;
  cols[ncols++] = globalRow(gi, gj, k);                          // center
  if (gi - 1 >= 0) cols[ncols++] = globalRow(gi - 1, gj, k);     // W
  if (gi + 1 < nx_) cols[ncols++] = globalRow(gi + 1, gj, k);    // E
  if (gj - 1 >= 0) cols[ncols++] = globalRow(gi, gj - 1, k);     // S
  if (gj + 1 < ny_) cols[ncols++] = globalRow(gi, gj + 1, k);    // N
  if (k - 1 >= 0) cols[ncols++] = globalRow(gi, gj, k - 1);      // Up (toward surface)
  if (k + 1 < nz_) cols[ncols++] = globalRow(gi, gj, k + 1);     // Down
}

std::pair<int, int> Decomp3D::ownershipRange() const {
  const int start = rowStartOfRank(comm_.px(), comm_.py());
  return {start, start + comm_.localNx() * comm_.localNy() * nz_};
}

void Decomp3D::ownedRowStencil(int local_owned, int& global_row, int* cols,
                               int& ncols) const {
  const int lnx = comm_.localNx();
  const int k = local_owned % nz_;
  const int horiz = local_owned / nz_;
  const int li = horiz % lnx;
  const int lj = horiz / lnx;
  const int gi = comm_.i0() + li;
  const int gj = comm_.j0() + lj;
  global_row = globalRow(gi, gj, k);
  stencilColumns(gi, gj, k, cols, ncols);
}

void Decomp3D::haloExchange(const RealArr1D& field) const {
  comm_.haloExchange3D(field, comm_.localNx() + 2, comm_.localNy() + 2, nz_);
}

void Decomp3D::packOwned(const RealArr1D& field_halo, const RealArr1D& owned) const {
  const int lnx = comm_.localNx();
  const int lny = comm_.localNy();
  const int nxh = lnx + 2;
  const int plane = nxh * (lny + 2);
  auto fh = Kokkos::create_mirror_view(field_halo);
  Kokkos::deep_copy(fh, field_halo);
  auto oh = Kokkos::create_mirror_view(owned);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      for (int k = 0; k < nz_; ++k)
        oh((lj * lnx + li) * nz_ + k) =
            fh((li + 1) + (lj + 1) * nxh + (k + 1) * plane);
  Kokkos::deep_copy(owned, oh);
}

void Decomp3D::unpackOwned(const RealArr1D& owned, const RealArr1D& field_halo) const {
  const int lnx = comm_.localNx();
  const int lny = comm_.localNy();
  const int nxh = lnx + 2;
  const int plane = nxh * (lny + 2);
  auto oh = Kokkos::create_mirror_view(owned);
  Kokkos::deep_copy(oh, owned);
  auto fh = Kokkos::create_mirror_view(field_halo);
  Kokkos::deep_copy(fh, field_halo);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      for (int k = 0; k < nz_; ++k)
        fh((li + 1) + (lj + 1) * nxh + (k + 1) * plane) =
            oh((lj * lnx + li) * nz_ + k);
  Kokkos::deep_copy(field_halo, fh);
}

}  // namespace frehg2
