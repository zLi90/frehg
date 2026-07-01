// P13.3.4 acceptance: non-uniform soil is wired through the production Orchestrator path
// (buildGroundwater -> buildSoilMap -> ReSolver::setSoilMap), reachable from run()/main().
//
// A GW-only 2-column box (no-flow, gravity redistribution from a uniform IC; lateral K=0 so
// columns are independent) is run two ways:
//   (A) uniform soil (no soil.map)              -> both columns use class 0.
//   (B) soil.map.from_file with class 0 | 1     -> column 0 = class 0, column 1 = class 1.
// Assertions:
//   - (B) column 0 reproduces (A) column 0 (< 1e-6): the class-0 path is unchanged.
//   - (B) column 0 and column 1 differ by a wide margin (> 1e-4): the per-cell map is actually
//     applied through the Orchestrator (a dropped/ignored map would make the columns identical).
#include <petscksp.h>

#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <mpi.h>
#include <string>
#include <vector>

#include "core/Grid.hpp"
#include "driver/Orchestrator.hpp"
#include "frehg2_test.hpp"
#include "io/Config.hpp"
#include "re/ReSolver.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;
constexpr int kNz = 20;

std::string gwSoilConfig(const std::string& out_h5, bool with_map, const std::string& map_file) {
  std::string s =
      "schema_version: '2.0'\n"
      "simulation: {id: gw_soil, mode: groundwater}\n"
      "domain: {nx: 2, ny: 1, nz: " + std::to_string(kNz) +
      ", dx: 1.0, dy: 1.0, dz: 0.01, dz_incre: 1.0, botz: -0.2, follow_terrain: false}\n"
      "time: {dt: 0.001, t_end: 1.0e9, max_steps: 100, output_interval: 0}\n"
      "modules: {surface_water: false, groundwater: true, solute: false}\n"
      "groundwater: {solver: pca, full_3d: false, adaptive_dt: true, use_corrector: true,\n"
      "  use_vg: true, use_mvg: false, air_entry_value: -0.02, dt_max: 2.0, dt_min: 0.001,\n"
      "  co_max: 2.0, specific_storage: 1.0e-5, bc_type_gw: [0,0,0,0,0,0],\n"
      "  solver: {ksp_type: cg, pc_type: jacobi, rtol: 1.0e-13}}\n"
      "soil:\n";
  if (with_map)
    s += "  map: {from_file: true, file: " + map_file + ", format: list}\n";
  s +=
      "  types:\n"
      "    - {id: 0, theta_s: 0.33, theta_r: 0.0, vg: {alpha: 1.43, n: 1.56}, k_sat: {x: 0.0, y: 0.0, z: 1.0e-6}}\n"
      "    - {id: 1, theta_s: 0.33, theta_r: 0.0, vg: {alpha: 1.43, n: 1.56}, k_sat: {x: 0.0, y: 0.0, z: 1.0e-4}}\n"
      "initial_conditions: {groundwater: {wc: 0.20}}\n"
      "output: {format: hdf5, filename: " + out_h5 + ", io_mode: serial_gather}\n";
  return s;
}

// Column water-content profile (k = 0..nz-1) for owned column i.
std::vector<double> columnProfile(const Orchestrator& orch, int i) {
  const Grid& g = orch.re()->grid();
  const auto& f = orch.re()->fields();
  std::vector<double> p(static_cast<size_t>(kNz));
  for (int k = 0; k < kNz; ++k) p[static_cast<size_t>(k)] = f.wc(g.getIndex(i, 0, k));
  return p;
}

double worstDiff(const std::vector<double>& a, const std::vector<double>& b) {
  double w = 0.0;
  for (size_t i = 0; i < a.size(); ++i) w = std::max(w, std::fabs(a[i] - b[i]));
  return w;
}
}  // namespace

TEST_CASE("orchestrator: non-uniform soil map is parsed, sliced, and applied via run()") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 1) return;

  // Write a 2-cell class-index list file (column 0 -> class 0, column 1 -> class 1).
  const std::string map_file = std::string(kTmp) + "/soil_class.txt";
  {
    std::ofstream o(map_file);
    o << "0\n1\n";
  }

  std::vector<double> uni_c0, map_c0, map_c1;
  {
    Config cfg = Config::fromString(
        gwSoilConfig(std::string(kTmp) + "/gw_soil_uniform.h5", false, ""), "");
    Orchestrator orch;
    orch.initialize(cfg);
    orch.run();
    uni_c0 = columnProfile(orch, 0);
  }
  {
    Config cfg = Config::fromString(
        gwSoilConfig(std::string(kTmp) + "/gw_soil_map.h5", true, map_file), "");
    Orchestrator orch;
    orch.initialize(cfg);
    orch.run();
    map_c0 = columnProfile(orch, 0);
    map_c1 = columnProfile(orch, 1);
  }

  // (1) class-0 column of the mapped run matches the uniform run's column 0.
  const double d0 = worstDiff(map_c0, uni_c0);
  std::fprintf(stderr, "  orch soil: worst |map_col0 - uniform_col0| = %.3e\n", d0);
  REQUIRE(d0 < 1.0e-6);

  // (2) class-0 and class-1 columns of the mapped run genuinely differ.
  const double d01 = worstDiff(map_c0, map_c1);
  std::fprintf(stderr, "  orch soil: worst |map_col0 - map_col1| = %.3e\n", d01);
  REQUIRE(d01 > 1.0e-4);
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
