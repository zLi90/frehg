// P2.4 acceptance: Decomp3D ownership, halo (z on-rank), 7-point stencil.
// Runs on 1 and 2 MPI ranks.
#include <mpi.h>

#include <vector>

#include <Kokkos_Core.hpp>

#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"

using namespace frehg2;

namespace {
void procGrid(int size, int& mpi_nx, int& mpi_ny) {
  if (size == 1) { mpi_nx = 1; mpi_ny = 1; }
  else if (size == 2) { mpi_nx = 2; mpi_ny = 1; }
  else if (size == 4) { mpi_nx = 2; mpi_ny = 2; }
  else { mpi_nx = size; mpi_ny = 1; }
}
}  // namespace

TEST_CASE("single-rank 3x2x4 ownership") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;
  MpiComm mc(3, 2, 1, 1);
  Decomp3D dd(mc, 4);
  auto [nx, ny, nz] = dd.globalSize();
  REQUIRE(nx == 3);
  REQUIRE(ny == 2);
  REQUIRE(nz == 4);
  REQUIRE(dd.globalRowCount() == 24);
  REQUIRE(dd.ownershipRange().first == 0);
  REQUIRE(dd.ownershipRange().second == 24);
}

TEST_CASE("globalRow bijection onto [0, nx*ny*nz)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  const int NX = 4, NY = 2, NZ = 3;
  MpiComm mc(NX, NY, pnx, pny);
  Decomp3D dd(mc, NZ);
  std::vector<int> hits(NX * NY * NZ, 0);
  for (int gj = 0; gj < NY; ++gj)
    for (int gi = 0; gi < NX; ++gi)
      for (int k = 0; k < NZ; ++k) {
        int r = dd.globalRow(gi, gj, k);
        REQUIRE(r >= 0);
        REQUIRE(r < NX * NY * NZ);
        hits[r]++;
      }
  for (int r = 0; r < NX * NY * NZ; ++r) REQUIRE(hits[r] == 1);
}

TEST_CASE("ownership ranges partition the global rows") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  const int NX = 4, NY = 2, NZ = 3;
  MpiComm mc(NX, NY, pnx, pny);
  Decomp3D dd(mc, NZ);
  auto [lo, hi] = dd.ownershipRange();
  std::vector<int> owned(NX * NY * NZ, 0);
  for (int r = lo; r < hi; ++r) owned[r] = 1;
  MPI_Allreduce(MPI_IN_PLACE, owned.data(), NX * NY * NZ, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  for (int r = 0; r < NX * NY * NZ; ++r) REQUIRE(owned[r] == 1);
}

TEST_CASE("7-point stencil: interior 7, vertical/edge fewer") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;
  MpiComm mc(3, 3, 1, 1);
  Decomp3D dd(mc, 3);
  int cols[kMaxStencil];
  int ncols = 0;
  dd.stencilColumns(1, 1, 1, cols, ncols);  // fully interior
  REQUIRE(ncols == 7);
  dd.stencilColumns(1, 1, 0, cols, ncols);  // top layer: no Up
  REQUIRE(ncols == 6);
  dd.stencilColumns(0, 0, 0, cols, ncols);  // corner top: W,S,Up missing
  REQUIRE(ncols == 4);
}

TEST_CASE("z neighbors are addressable on-rank (no MPI in z)") {
  // The vertical column neighbors of an owned cell share the same rank: globalRow(.,.,k)
  // and globalRow(.,.,k+1) are both in this rank's ownership range.
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(4, 2, pnx, pny);
  Decomp3D dd(mc, 5);
  auto [lo, hi] = dd.ownershipRange();
  const int gi = mc.i0();
  const int gj = mc.j0();
  for (int k = 0; k < 4; ++k) {
    int r0 = dd.globalRow(gi, gj, k);
    int r1 = dd.globalRow(gi, gj, k + 1);
    REQUIRE(r0 >= lo);
    REQUIRE(r0 < hi);
    REQUIRE(r1 >= lo);
    REQUIRE(r1 < hi);
  }
}

TEST_CASE("pack/unpack owned round-trip (3D)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(4, 2, pnx, pny);
  const int NZ = 5;
  Decomp3D dd(mc, NZ);

  RealArr1D field("field", dd.localHaloFieldSize());
  RealArr1D owned("owned", dd.ownedRowCount());
  auto fh = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(fh, 0.0);
  const int lnx = mc.localNx(), lny = mc.localNy(), nxh = lnx + 2;
  const int plane = nxh * (lny + 2);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      for (int k = 0; k < NZ; ++k)
        fh((li + 1) + (lj + 1) * nxh + (k + 1) * plane) = static_cast<real>(1000 * k + 10 * lj + li);
  Kokkos::deep_copy(field, fh);

  dd.packOwned(field, owned);
  RealArr1D field2("field2", dd.localHaloFieldSize());
  dd.unpackOwned(owned, field2);
  auto f2 = Kokkos::create_mirror_view(field2);
  Kokkos::deep_copy(f2, field2);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      for (int k = 0; k < NZ; ++k)
        REQUIRE(f2((li + 1) + (lj + 1) * nxh + (k + 1) * plane) ==
                Approx(1000 * k + 10 * lj + li).margin(1e-12));
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int rc = frehg2test::runAll();
  Kokkos::finalize();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return grc;
}
