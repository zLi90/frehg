// P4.6 (4c) b1-sw benchmark gate, serial one rank.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "frehg2_test.hpp"
#include "swe/b1_sw_runner.hpp"

using namespace frehg2;
using namespace frehg2::b1_sw;

TEST_CASE("b1-sw serial: per-cell depth matches legacy reference (relative L2 < 1e-3)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  B1Params params = loadB1Params(FREHG2_LEGACY_B1_DIR);
  MpiComm mc(params.gnx, params.gny, 1, 1);
  B1RunResult result = runB1(params, mc, FREHG2_LEGACY_B1_DIR);

  std::fprintf(stderr,
               "  b1-sw aggregate time-series relative L2 = %.3e "
               "(worst single snapshot = %.3e at t=%d)\n",
               result.aggregate_rel_l2_vs_legacy, result.worst_snapshot_rel_l2,
               result.worst_t);
  REQUIRE(result.aggregate_rel_l2_vs_legacy < 1.0e-3);
  REQUIRE(result.worst_snapshot_rel_l2 < 5.0e-3);
}

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  MPI_Init(&argc, &argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  int rc = frehg2test::runAll();
  int grc = rc;
  MPI_Allreduce(MPI_IN_PLACE, &grc, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  PetscFinalize();
  MPI_Finalize();
  Kokkos::finalize();
  return grc;
}
