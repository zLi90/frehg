// P10 Task 10.3.5 — b1-sw end-to-end on the GPU backend.
//
// On a Kokkos-CUDA/HIP build (FREHG2_GPU_ASSEMBLY==1), this re-runs the blocking b1-sw gate
// through the GPU-native assembly path (MatSetValuesCOO + device VECKOKKOS bridge) and asserts
// the solution still matches the legacy reference — i.e. the GPU result equals the validated
// CPU/OpenMP result within the gate tolerance. The b1-sw runner is the SAME one the CPU gate
// uses, so the only thing that changes here is the execution space and the device bridge.
//
// This test is labeled `gpu` and is DISABLED on macOS (no CUDA): the body below compiles to a
// skip, and CTest never runs it locally (see tests/gpu/CMakeLists.txt). Real execution is the
// deferred Linux/NVIDIA validation per docs/gpu_validation_policy.md.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "frehg2_test.hpp"
#include "linear/backends/GpuAssembly.hpp"
#include "swe/b1_sw_runner.hpp"

using namespace frehg2;

TEST_CASE("b1-sw on GPU matches the legacy reference (deferred Linux/NVIDIA)") {
#if FREHG2_GPU_ASSEMBLY
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;
  b1_sw::B1Params params = b1_sw::loadB1Params(FREHG2_LEGACY_B1_DIR);
  MpiComm mc(params.gnx, params.gny, 1, 1);
  b1_sw::B1RunResult result = b1_sw::runB1(params, mc, FREHG2_LEGACY_B1_DIR);
  std::fprintf(stderr, "  [gpu] b1-sw aggregate rel-L2 vs legacy = %.3e\n",
               result.aggregate_rel_l2_vs_legacy);
  // Same gate as the CPU path (user-authorized 1e-3); GPU must reproduce the CPU/OpenMP result.
  REQUIRE(result.aggregate_rel_l2_vs_legacy < 1.0e-3);
  REQUIRE(result.worst_snapshot_rel_l2 < 5.0e-3);
#else
  // No GPU backend in this build: skipped. The CPU gate (tests/swe/test_swe_b1.cpp) is the
  // active b1-sw validation here.
  REQUIRE(true);
#endif
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
