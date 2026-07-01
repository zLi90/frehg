// P2.4 acceptance: DomainDecomposition is backend-agnostic. This translation unit and its
// target link ONLY Kokkos + MPI (frehg2_decomp), with NO PETSc. If DomainDecomposition
// pulled in any PETSc symbol, this target would fail to link. It exercises the public API
// to ensure it is fully usable without a solver library.
#include <mpi.h>

#include <Kokkos_Core.hpp>

#include "frehg2_test.hpp"
#include "linear/DomainDecomposition.hpp"

using namespace frehg2;

TEST_CASE("Decomp2D usable with no PETSc linked") {
  MpiComm mc(5, 4, 1, 1);
  Decomp2D dd(mc);
  REQUIRE(dd.globalRowCount() == 20);
  int cols[kMaxStencil];
  int ncols = 0;
  dd.stencilColumns(2, 2, cols, ncols);
  REQUIRE(ncols == 5);

  RealArr1D field("f", dd.localHaloFieldSize());
  RealArr1D owned("o", dd.ownedRowCount());
  Kokkos::deep_copy(field, 1.0);
  dd.packOwned(field, owned);
  auto oh = Kokkos::create_mirror_view(owned);
  Kokkos::deep_copy(oh, owned);
  REQUIRE(oh(0) == Approx(1.0).margin(1e-12));
}

TEST_CASE("Decomp3D usable with no PETSc linked") {
  MpiComm mc(3, 3, 1, 1);
  Decomp3D dd(mc, 4);
  REQUIRE(dd.globalRowCount() == 36);
  int cols[kMaxStencil];
  int ncols = 0;
  dd.stencilColumns(1, 1, 1, cols, ncols);
  REQUIRE(ncols == 7);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int rc = frehg2test::runAll();
  Kokkos::finalize();
  MPI_Finalize();
  return rc;
}
