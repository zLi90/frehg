// P10 Task 10.3.7 — GPU variant of the async SW<->GW coupling equivalence gate.
//
// On a Kokkos-CUDA/HIP build (FREHG2_GPU_ASSEMBLY==1), this runs a coupled scenario in both
// coupling.mode=sync and coupling.mode=async on the GPU backend and asserts they agree to
// 1e-10 (the P11 gate carried forward). The async GPU realization (PetscSubcomm split +
// CUDA/HIP streams, src/linear/backends/PetscSubcommSplit) keeps the SAME Gauss-Seidel
// operation sequence and exchange point as the synchronous path, so the result is unchanged;
// see docs/research_notes/async_gpu.md.
//
// Labeled `gpu`, DISABLED on macOS (no CUDA): the body compiles to a skip and never runs
// locally. The CPU async equivalence (tests/integration/test_async_coupling.cpp) is the active
// numerical oracle here. Deferred to Linux/NVIDIA per docs/gpu_validation_policy.md.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <string>

#include "frehg2/core/define.hpp"
#include "frehg2_test.hpp"
#include "linear/backends/GpuAssembly.hpp"

#if FREHG2_GPU_ASSEMBLY
#include "driver/Orchestrator.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "re/ReSolver.hpp"
#include "swe/SweSolver.hpp"

namespace {
using namespace frehg2;
RealArr1DHost ownedWc(const Orchestrator& o) {
  const ReSolver* re = o.re();
  const auto& f = re->fields();
  const int nx = re->grid().nx(), ny = re->grid().ny(), nz = re->grid().nz();
  RealArr1DHost out("wc", static_cast<size_t>(nx * ny * nz));
  size_t e = 0;
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) out(e++) = f.wc(re->grid().getIndex(i, j, k));
  return out;
}
}  // namespace
#endif

using namespace frehg2;

TEST_CASE("coupled sync vs async on GPU agree to 1e-10 (deferred Linux/NVIDIA)") {
#if FREHG2_GPU_ASSEMBLY
  using namespace frehg2::orch_test;
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;  // async GPU path is single-rank per window; multi-rank uses subcomm

  const int nx = 4, ny = 4, nz = 6;
  const double dt = 1.0, t_end = 20.0;
  const long long nsteps = 20;
  const double init_eta = 0.1, init_wc = 0.2;
  const std::string base = std::string(FREHG2_IO_TMP);

  Config cs = Config::fromString(
      coupledConfig(base + "/gpu_sync.h5", nx, ny, nz, dt, t_end, nsteps, init_eta, init_wc,
                    "sync"),
      "");
  Orchestrator osync;
  osync.initialize(cs);
  osync.run();
  RealArr1DHost wc_sync = ownedWc(osync);

  Config ca = Config::fromString(
      coupledConfig(base + "/gpu_async.h5", nx, ny, nz, dt, t_end, nsteps, init_eta, init_wc,
                    "async"),
      "");
  Orchestrator oasync;
  oasync.initialize(ca);
  oasync.run();
  RealArr1DHost wc_async = ownedWc(oasync);

  const double l2 = relL2(wc_async, wc_sync);
  std::fprintf(stderr, "  [gpu] coupled async-vs-sync rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-10);
#else
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
