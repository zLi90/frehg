// P10 Task 10.3.5 — b2-gw end-to-end on the GPU backend.
//
// On a Kokkos-CUDA/HIP build (FREHG2_GPU_ASSEMBLY==1), this runs the Richards-equation b2-gw
// column twice on the GPU — once through the production Orchestrator and once through the direct
// RE path (exactly as tests/integration/test_orchestrator_parity.cpp does on CPU) — and asserts
// they agree to L2 < 1e-6 while going through the GPU-native COO matrix assembly and device Vec
// bridge. Combined with test_b1_sw_gpu (which compares the GPU result to the legacy reference),
// this validates the GW solve on the device path.
//
// Labeled `gpu`, DISABLED on macOS (no CUDA): the body compiles to a skip and never runs
// locally. Deferred to Linux/NVIDIA per docs/gpu_validation_policy.md.
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <string>

#include "frehg2/core/define.hpp"
#include "frehg2_test.hpp"
#include "linear/backends/GpuAssembly.hpp"

#if FREHG2_GPU_ASSEMBLY
#include "core/Grid.hpp"
#include "core/MpiComm.hpp"
#include "driver/Orchestrator.hpp"
#include "integration/orch_test_util.hpp"
#include "io/Config.hpp"
#include "linear/DomainDecomposition.hpp"
#include "linear/backends/PetscLinearSolver.hpp"
#include "re/ReSolver.hpp"
#endif

using namespace frehg2;

TEST_CASE("b2-gw on GPU: Orchestrator reproduces the direct RE path (deferred Linux/NVIDIA)") {
#if FREHG2_GPU_ASSEMBLY
  using namespace frehg2::orch_test;
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  const int nz = 50;
  const double dt = 1.0e-4;
  const long long nsteps = 50;
  const double t_end = 1.0e6;
  const double init_wc = 0.033;

  MpiComm mc(1, 1, 1, 1);
  Grid grid(1, 1, nz, 1.0, 1.0, 0.01, 1.0);
  ReSolver re(grid, &mc);
  ReParams rp;
  rp.soil.alpha = 1.43;
  rp.soil.n = 1.56;
  rp.soil.theta_s = 0.33;
  rp.soil.theta_r = 0.0;
  rp.soil.Ks_z = 2.89e-6;
  rp.soil.Ss = 1.0e-5;
  rp.soil.use_vg = true;
  rp.dx = 1.0;
  rp.dy = 1.0;
  rp.dz = 0.01;
  rp.botz = -1.0;
  rp.dt = dt;
  rp.dt_min = dt;
  rp.dt_max = 2.0;
  rp.co_max = 2.0;
  rp.adaptive_dt = true;
  rp.use_corrector = true;
  rp.bc_type[5] = 1;
  rp.htop = 0.0;
  re.setParams(rp);
  re.initializeUniformColumn(init_wc);
  Decomp3D dd(mc, nz);
  SolverConfig cfg;
  cfg.ksp_type = "cg";
  cfg.pc_type = "jacobi";
  cfg.rtol = 1.0e-10;
  PetscLinearSolver solver(cfg);
  re.attachSolver(solver, dd);
  for (long long step = 0; step < nsteps; ++step) re.advanceStep();
  RealArr1DHost direct("direct", static_cast<size_t>(nz));
  for (int k = 0; k < nz; ++k) direct(static_cast<size_t>(k)) = re.fields().wc(grid.getIndex(0, 0, k));

  const std::string out = std::string(FREHG2_IO_TMP) + "/gpu_b2.h5";
  Config c = Config::fromString(gwConfig(out, nz, dt, t_end, nsteps, init_wc), "");
  Orchestrator orch;
  orch.initialize(c);
  orch.run();
  RealArr1DHost driver("driver", static_cast<size_t>(nz));
  for (int k = 0; k < nz; ++k)
    driver(static_cast<size_t>(k)) = orch.re()->fields().wc(grid.getIndex(0, 0, k));

  const double l2 = relL2(driver, direct);
  std::fprintf(stderr, "  [gpu] b2-gw orchestrator-vs-direct rel-L2 = %.3e\n", l2);
  REQUIRE(l2 < 1.0e-6);
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
