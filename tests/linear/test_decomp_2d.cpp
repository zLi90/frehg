// P2.4 acceptance: Decomp2D ownership, bijective global numbering, halo, stencil.
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

TEST_CASE("single-rank 4x3 ownership and corners") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // this case is 1-rank only
  MpiComm mc(4, 3, 1, 1);
  Decomp2D dd(mc);
  REQUIRE(dd.globalSize().first == 4);
  REQUIRE(dd.globalSize().second == 3);
  REQUIRE(dd.ownershipRange().first == 0);
  REQUIRE(dd.ownershipRange().second == 12);
  auto [i0, j0, ni, nj] = dd.localCorners();
  REQUIRE(i0 == 0);
  REQUIRE(j0 == 0);
  REQUIRE(ni == 4);
  REQUIRE(nj == 3);
}

TEST_CASE("globalRow is a bijection onto [0, nx*ny)") {
  int size = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  const int NX = 6, NY = 4;
  MpiComm mc(NX, NY, pnx, pny);
  Decomp2D dd(mc);

  // Every rank computes globalRow for all cells; the union must hit each row once.
  std::vector<int> hits(NX * NY, 0);
  for (int gj = 0; gj < NY; ++gj)
    for (int gi = 0; gi < NX; ++gi) {
      int r = dd.globalRow(gi, gj);
      REQUIRE(r >= 0);
      REQUIRE(r < NX * NY);
      hits[r]++;
    }
  // Each rank sees the full domain identically, so each row is hit exactly once here.
  for (int r = 0; r < NX * NY; ++r) REQUIRE(hits[r] == 1);
}

TEST_CASE("ownership ranges partition [0, nx*ny) with no gaps/overlaps") {
  int size = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  const int NX = 6, NY = 4;
  MpiComm mc(NX, NY, pnx, pny);
  Decomp2D dd(mc);

  auto [lo, hi] = dd.ownershipRange();
  std::vector<int> owned_count(NX * NY, 0);
  for (int r = lo; r < hi; ++r) owned_count[r] = 1;
  MPI_Allreduce(MPI_IN_PLACE, owned_count.data(), NX * NY, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  for (int r = 0; r < NX * NY; ++r) REQUIRE(owned_count[r] == 1);
}

TEST_CASE("stencilColumns: interior 5, edge fewer") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;
  MpiComm mc(5, 5, 1, 1);
  Decomp2D dd(mc);
  int cols[kMaxStencil];
  int ncols = 0;
  dd.stencilColumns(2, 2, cols, ncols);
  REQUIRE(ncols == 5);
  dd.stencilColumns(0, 0, cols, ncols);  // corner
  REQUIRE(ncols == 3);
  dd.stencilColumns(0, 2, cols, ncols);  // left edge
  REQUIRE(ncols == 4);
}

TEST_CASE("haloExchange fills ghosts with neighbor-owned global rows' values") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  const int NX = 8, NY = 6;
  MpiComm mc(NX, NY, pnx, pny);
  Decomp2D dd(mc);

  const int lnx = mc.localNx();
  const int lny = mc.localNy();
  const int nxh = lnx + 2;
  RealArr1D field("field", dd.localHaloFieldSize());
  auto fh = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(fh, -1.0);
  // store the global row index as the value in each owned cell
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      fh((li + 1) + (lj + 1) * nxh) =
          static_cast<real>(dd.globalRow(mc.i0() + li, mc.j0() + lj));
  Kokkos::deep_copy(field, fh);

  dd.haloExchange(field);
  Kokkos::deep_copy(fh, field);

  // left ghost should equal the global row of the neighbor cell (gi-1,gj)
  if (mc.leftRank() != MPI_PROC_NULL) {
    for (int lj = 0; lj < lny; ++lj)
      REQUIRE(fh(0 + (lj + 1) * nxh) ==
              Approx(static_cast<real>(dd.globalRow(mc.i0() - 1, mc.j0() + lj)))
                  .margin(1e-9));
  }
}

TEST_CASE("pack/unpack owned round-trip") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(8, 6, pnx, pny);
  Decomp2D dd(mc);

  RealArr1D field("field", dd.localHaloFieldSize());
  RealArr1D owned("owned", dd.ownedRowCount());
  auto fh = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(fh, 0.0);
  const int lnx = mc.localNx(), lny = mc.localNy(), nxh = lnx + 2;
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      fh((li + 1) + (lj + 1) * nxh) = static_cast<real>(100 * lj + li);
  Kokkos::deep_copy(field, fh);

  dd.packOwned(field, owned);
  RealArr1D field2("field2", dd.localHaloFieldSize());
  dd.unpackOwned(owned, field2);

  auto f2 = Kokkos::create_mirror_view(field2);
  Kokkos::deep_copy(f2, field2);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      REQUIRE(f2((li + 1) + (lj + 1) * nxh) == Approx(100 * lj + li).margin(1e-12));
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
