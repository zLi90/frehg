// P2.3 acceptance: MpiComm decomposition queries, halo exchange, gather.
// Runs on 1 and 4 MPI ranks (custom main initializes MPI + Kokkos).
#include <mpi.h>

#include <Kokkos_Core.hpp>

#include "core/MpiComm.hpp"
#include "frehg2_test.hpp"

using namespace frehg2;

namespace {
// Distinct value per global cell.
real fval(int gi, int gj) { return static_cast<real>(gi) + 1000.0 * static_cast<real>(gj); }

// Pick a process grid for the current comm size.
void procGrid(int size, int& mpi_nx, int& mpi_ny) {
  if (size == 1) { mpi_nx = 1; mpi_ny = 1; }
  else if (size == 2) { mpi_nx = 2; mpi_ny = 1; }
  else if (size == 4) { mpi_nx = 2; mpi_ny = 2; }
  else { mpi_nx = size; mpi_ny = 1; }
}
}  // namespace

TEST_CASE("ownership and global<->local round-trip") {
  int size = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(10, 10, pnx, pny);

  for (int lj = 0; lj < mc.localNy(); ++lj) {
    for (int li = 0; li < mc.localNx(); ++li) {
      auto [gi, gj] = mc.localToGlobal(li, lj);
      auto [li2, lj2, owner] = mc.globalToLocal(gi, gj);
      REQUIRE(li2 == li);
      REQUIRE(lj2 == lj);
      REQUIRE(owner == rank);
    }
  }
}

TEST_CASE("halo exchange fills ghosts from neighbor interiors") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  MpiComm mc(10, 10, pnx, pny);

  const int lnx = mc.localNx();
  const int lny = mc.localNy();
  const int nxh = lnx + 2;
  const int nyh = lny + 2;

  RealArr1D field("field", static_cast<size_t>(nxh) * nyh);
  auto h = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(h, 0.0);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      h((li + 1) + (lj + 1) * nxh) = fval(mc.i0() + li, mc.j0() + lj);
  Kokkos::deep_copy(field, h);

  mc.haloExchange2D(field, nxh, nyh);

  Kokkos::deep_copy(h, field);
  // left ghost
  if (mc.leftRank() != MPI_PROC_NULL) {
    for (int lj = 0; lj < lny; ++lj)
      REQUIRE(h(0 + (lj + 1) * nxh) == Approx(fval(mc.i0() - 1, mc.j0() + lj)).margin(1e-12));
  }
  if (mc.rightRank() != MPI_PROC_NULL) {
    for (int lj = 0; lj < lny; ++lj)
      REQUIRE(h((lnx + 1) + (lj + 1) * nxh) ==
              Approx(fval(mc.i0() + lnx, mc.j0() + lj)).margin(1e-12));
  }
  if (mc.downRank() != MPI_PROC_NULL) {
    for (int li = 0; li < lnx; ++li)
      REQUIRE(h((li + 1) + 0 * nxh) == Approx(fval(mc.i0() + li, mc.j0() - 1)).margin(1e-12));
  }
  if (mc.upRank() != MPI_PROC_NULL) {
    for (int li = 0; li < lnx; ++li)
      REQUIRE(h((li + 1) + (lny + 1) * nxh) ==
              Approx(fval(mc.i0() + li, mc.j0() + lny)).margin(1e-12));
  }
}

TEST_CASE("gatherToRank0 reassembles the full grid") {
  int size = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int pnx = 1, pny = 1;
  procGrid(size, pnx, pny);
  const int NX = 10, NY = 10;
  MpiComm mc(NX, NY, pnx, pny);

  const int lnx = mc.localNx();
  const int lny = mc.localNy();
  const int nxh = lnx + 2;
  const int nyh = lny + 2;
  RealArr1D field("field", static_cast<size_t>(nxh) * nyh);
  auto h = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(h, 0.0);
  for (int lj = 0; lj < lny; ++lj)
    for (int li = 0; li < lnx; ++li)
      h((li + 1) + (lj + 1) * nxh) = fval(mc.i0() + li, mc.j0() + lj);
  Kokkos::deep_copy(field, h);

  RealArr1DHost global("global", static_cast<size_t>(NX) * NY);
  Kokkos::deep_copy(global, -1.0);
  mc.gatherToRank0(field, nxh, nyh, global);

  if (rank == 0) {
    for (int gj = 0; gj < NY; ++gj)
      for (int gi = 0; gi < NX; ++gi)
        REQUIRE(global(gi + gj * NX) == Approx(fval(gi, gj)).margin(1e-12));
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int rc = frehg2test::runAll();
  Kokkos::finalize();
  int global_rc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &global_rc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return global_rc;
}
