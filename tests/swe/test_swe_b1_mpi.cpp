// P4.4e: b1-sw MPI rank-count equivalence gate.
//
// Runs the full legacy-exact b1-sw case on 2 and 4 MPI ranks (y-decomposition of the
// 1x10 channel), gathers global depth fields, and compares against the embedded 1-rank
// Frehg2 reference snapshots (deterministic serial run). Gate: relative L2 < 1e-10.
#include <petscksp.h>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include "frehg2_test.hpp"
#include "swe/b1_sw_runner.hpp"

using namespace frehg2;
using namespace frehg2::b1_sw;

namespace {
const char* kLegacy = FREHG2_LEGACY_B1_DIR;
const char* kNp1Ref = FREHG2_B1_NP1_REF;

std::vector<B1Snapshot> loadNp1Reference() {
  std::ifstream in(kNp1Ref);
  if (!in) throw std::runtime_error("cannot open np1 reference: " + std::string(kNp1Ref));
  std::vector<B1Snapshot> snaps;
  int t = 0;
  while (in >> t) {
    B1Snapshot s;
    s.t = t;
    s.depth.resize(10);
    for (int j = 0; j < 10; ++j) in >> s.depth[static_cast<size_t>(j)];
    snaps.push_back(std::move(s));
  }
  return snaps;
}

double maxRelL2VsReference(const std::vector<B1Snapshot>& got,
                           const std::vector<B1Snapshot>& ref) {
  REQUIRE(got.size() == ref.size());
  double worst = 0.0;
  for (size_t k = 0; k < got.size(); ++k) {
    REQUIRE(got[k].t == ref[k].t);
    REQUIRE(got[k].depth.size() == ref[k].depth.size());
    worst = std::max(worst, relL2(got[k].depth, ref[k].depth));
  }
  return worst;
}
}  // namespace

TEST_CASE("b1-sw MPI: 2-rank y-decomp matches 1-rank reference (L2 < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) return;

  int pnx = 1, pny = 2;
  B1Params params = loadB1Params(kLegacy);
  MpiComm mc(params.gnx, params.gny, pnx, pny);
  B1RunResult result = runB1(params, mc, kLegacy);

  if (mc.rank() == 0) {
    const auto ref = loadNp1Reference();
    const double worst = maxRelL2VsReference(result.snapshots, ref);
    std::fprintf(stderr, "  b1-sw np2 vs np1 reference: worst snapshot rel-L2 = %.3e\n",
                 worst);
    REQUIRE(worst < 1.0e-10);
  }
}

TEST_CASE("b1-sw MPI: 4-rank y-decomp matches 1-rank reference (L2 < 1e-10)") {
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 4) return;

  int pnx = 1, pny = 4;
  B1Params params = loadB1Params(kLegacy);
  MpiComm mc(params.gnx, params.gny, pnx, pny);
  B1RunResult result = runB1(params, mc, kLegacy);

  if (mc.rank() == 0) {
    const auto ref = loadNp1Reference();
    const double worst = maxRelL2VsReference(result.snapshots, ref);
    std::fprintf(stderr, "  b1-sw np4 vs np1 reference: worst snapshot rel-L2 = %.3e\n",
                 worst);
    REQUIRE(worst < 1.0e-10);
  }
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
