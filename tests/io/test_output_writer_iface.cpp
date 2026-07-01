// P3.2.0 acceptance: solver-facing code uses only the OutputWriter interface; swapping the
// concrete writer needs no change in that code. The factory returns an Hdf5Writer here.
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <string>

#include "core/Grid.hpp"
#include "frehg2/io/OutputWriter.hpp"
#include "frehg2_test.hpp"
#include "io/Hdf5Reader.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;
constexpr int NX = 3, NY = 2;

// "Solver-facing" routine: knows nothing about HDF5, only OutputWriter.
void emitResults(OutputWriter& w, const Grid& grid, const RealArr1DHost& surf) {
  RunMetadata meta;
  meta.title = "iface test";
  w.openFile(std::string(kTmp) + "/iface.h5", meta);
  w.writeDomain(grid);
  w.writeSurfaceField("water_depth", 0.0, surf, "");
  w.writeMonitorSeries("gauge0", {0.0, 1.0, 2.0}, {3.0, 4.0, 5.0});
  w.close();
}
}  // namespace

TEST_CASE("factory yields an OutputWriter; solver code stays backend-agnostic") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  IoLayout L;
  L.nx = NX;
  L.ny = NY;
  L.nz = 1;
  // Single-rank owns all surface cells in this iface smoke test.
  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int total = NX * NY;
  int base = total / size, rem = total % size;
  int count = base + (rank < rem ? 1 : 0);
  int start = rank * base + (rank < rem ? rank : rem);
  L.surf_global_idx = IntArr1DHost("idx", static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) L.surf_global_idx(static_cast<size_t>(i)) = start + i;
  L.subs_global_idx = IntArr1DHost("idx3", 0);

  RealArr1DHost surf("surf", static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) surf(static_cast<size_t>(i)) = 10.0 + (start + i);

  Grid grid(NX, NY, 1, 1.0, 1.0, 1.0);

  std::unique_ptr<OutputWriter> w =
      makeOutputWriter("hdf5", MPI_COMM_WORLD, L, IoMode::SerialGather);
  REQUIRE(w != nullptr);
  emitResults(*w, grid, surf);

  if (rank == 0) {
    auto s = Hdf5Reader::readSurface(std::string(kTmp) + "/iface.h5", "water_depth", 0);
    REQUIRE(s.extent(0) == static_cast<size_t>(NX * NY));
    for (int g = 0; g < NX * NY; ++g) {
      REQUIRE(s(static_cast<size_t>(g)) == Approx(10.0 + g).margin(1e-9));
    }
  }
}

TEST_CASE("unknown format throws") {
  IoLayout L;
  bool threw = false;
  try {
    (void)makeOutputWriter("netcdf", MPI_COMM_WORLD, L, IoMode::SerialGather);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  REQUIRE(threw);
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
