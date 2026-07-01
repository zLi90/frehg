// P3.2 acceptance: HDF5 writer/reader round-trip, provenance metadata, multi-timestep, all
// three io_modes (serial_gather / parallel_collective / file_per_rank), and XDMF sidecar.
// Custom main initializes MPI + Kokkos; runs on 1 and 2 ranks.
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <fstream>
#include <sstream>
#include <string>

#include "core/Grid.hpp"
#include "frehg2/io/OutputWriter.hpp"
#include "frehg2_test.hpp"
#include "io/Hdf5Reader.hpp"
#include "io/Hdf5Writer.hpp"
#include "io/Sha256.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;  // CMAKE_CURRENT_BINARY_DIR

constexpr int NX = 4, NY = 3, NZ = 2;

void partition(int n, int size, int rank, int& start, int& count) {
  const int base = n / size;
  const int rem = n % size;
  count = base + (rank < rem ? 1 : 0);
  start = rank * base + (rank < rem ? rank : rem);
}

// Build an IoLayout where each rank owns a contiguous chunk of the global physical indices.
IoLayout makeLayout(MPI_Comm comm) {
  int size = 1, rank = 0;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
  IoLayout L;
  L.nx = NX;
  L.ny = NY;
  L.nz = NZ;
  int s = 0, c = 0;
  partition(NX * NY, size, rank, s, c);
  L.surf_global_idx = IntArr1DHost("surf_idx", static_cast<size_t>(c));
  for (int i = 0; i < c; ++i) L.surf_global_idx(static_cast<size_t>(i)) = s + i;
  int s3 = 0, c3 = 0;
  partition(NX * NY * NZ, size, rank, s3, c3);
  L.subs_global_idx = IntArr1DHost("subs_idx", static_cast<size_t>(c3));
  for (int i = 0; i < c3; ++i) L.subs_global_idx(static_cast<size_t>(i)) = s3 + i;
  return L;
}

double fsurf(int g, double t) { return 1.0 + 0.5 * g + 0.001 * t; }
double fsub(int g) { return 100.0 - 0.25 * g; }

RealArr1DHost ownedSurf(const IoLayout& L, double t) {
  RealArr1DHost v("ov", L.surf_global_idx.extent(0));
  for (size_t i = 0; i < v.extent(0); ++i)
    v(i) = static_cast<real>(fsurf(L.surf_global_idx(i), t));
  return v;
}
RealArr1DHost ownedSub(const IoLayout& L) {
  RealArr1DHost v("ovs", L.subs_global_idx.extent(0));
  for (size_t i = 0; i < v.extent(0); ++i)
    v(i) = static_cast<real>(fsub(L.subs_global_idx(i)));
  return v;
}

RunMetadata makeMeta() {
  RunMetadata m;
  m.title = "test run";
  m.frehg2_version = "0.1.0";
  m.git_sha = "abc123def456";
  m.git_dirty = true;
  m.config_sha256 = sha256Hex("resolved-config-text");
  m.build_type = "Debug";
  m.compiler = "GNU 15";
  m.kokkos_backend = "Serial";
  m.solver_backend = "petsc";
  return m;
}

std::string slurp(const std::string& path) {
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}
}  // namespace

TEST_CASE("serial_gather: round-trip, provenance, multi-timestep, xdmf") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  IoLayout L = makeLayout(MPI_COMM_WORLD);
  const std::string path = std::string(kTmp) + "/out_serial.h5";
  Grid grid(NX, NY, NZ, 80.0, 80.0, 0.1);

  {
    Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::SerialGather);
    w.openFile(path, makeMeta());
    w.writeDomain(grid);
    w.writeSurfaceField("water_depth", 0.0, ownedSurf(L, 0.0), "");
    w.writeSurfaceField("water_depth", 1800.0, ownedSurf(L, 1800.0), "");
    w.writeSubsurfaceField("hydraulic_head", 0.0, ownedSub(L), "");
    w.close();
  }

  if (rank == 0) {
    auto s0 = Hdf5Reader::readSurface(path, "water_depth", 0);
    auto s1 = Hdf5Reader::readSurface(path, "water_depth", 1800);
    REQUIRE(s0.extent(0) == static_cast<size_t>(NX * NY));
    for (int g = 0; g < NX * NY; ++g) {
      REQUIRE(s0(static_cast<size_t>(g)) == Approx(fsurf(g, 0.0)).margin(1e-9));
      REQUIRE(s1(static_cast<size_t>(g)) == Approx(fsurf(g, 1800.0)).margin(1e-9));
    }
    auto sub = Hdf5Reader::readSubsurface(path, "hydraulic_head", 0);
    REQUIRE(sub.extent(0) == static_cast<size_t>(NX * NY * NZ));
    for (int g = 0; g < NX * NY * NZ; ++g) {
      REQUIRE(sub(static_cast<size_t>(g)) == Approx(fsub(g)).margin(1e-9));
    }

    // Provenance.
    REQUIRE(Hdf5Reader::readMetadataAttr(path, "git_sha") == std::string("abc123def456"));
    REQUIRE(Hdf5Reader::readMetadataAttr(path, "config_sha256") ==
            sha256Hex("resolved-config-text"));
    REQUIRE(Hdf5Reader::readMetadataAttr(path, "io_mode") == std::string("serial_gather"));
    REQUIRE(Hdf5Reader::readSurfaceUnits(path, "water_depth", 0) == std::string("m"));

    // XDMF sidecar smoke check: exists, valid-ish XML, references the dataset.
    const std::string xmf = slurp(std::string(kTmp) + "/out_serial.xmf");
    REQUIRE(xmf.find("<Xdmf") != std::string::npos);
    REQUIRE(xmf.find("/surface/water_depth/0") != std::string::npos);
    REQUIRE(xmf.find("/surface/water_depth/1800") != std::string::npos);
  }
}

TEST_CASE("parallel_collective: shared file reassembles to global field") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  IoLayout L = makeLayout(MPI_COMM_WORLD);
  const std::string path = std::string(kTmp) + "/out_parallel.h5";
  Grid grid(NX, NY, NZ, 80.0, 80.0, 0.1);

  {
    Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::ParallelCollective);
    w.openFile(path, makeMeta());
    w.writeDomain(grid);
    w.writeSurfaceField("water_depth", 0.0, ownedSurf(L, 0.0), "");
    w.writeSubsurfaceField("hydraulic_head", 0.0, ownedSub(L), "");
    w.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    auto s0 = Hdf5Reader::readSurface(path, "water_depth", 0);
    REQUIRE(s0.extent(0) == static_cast<size_t>(NX * NY));
    for (int g = 0; g < NX * NY; ++g) {
      REQUIRE(s0(static_cast<size_t>(g)) == Approx(fsurf(g, 0.0)).margin(1e-9));
    }
    auto sub = Hdf5Reader::readSubsurface(path, "hydraulic_head", 0);
    for (int g = 0; g < NX * NY * NZ; ++g) {
      REQUIRE(sub(static_cast<size_t>(g)) == Approx(fsub(g)).margin(1e-9));
    }
    REQUIRE(Hdf5Reader::readMetadataAttr(path, "io_mode") ==
            std::string("parallel_collective"));
  }
}

TEST_CASE("file_per_rank: per-rank files reassemble to global field") {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  IoLayout L = makeLayout(MPI_COMM_WORLD);
  const std::string stem = std::string(kTmp) + "/out_fpr";
  const std::string path = stem + ".h5";
  Grid grid(NX, NY, NZ, 80.0, 80.0, 0.1);

  {
    Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::FilePerRank);
    w.openFile(path, makeMeta());
    w.writeDomain(grid);
    w.writeSurfaceField("water_depth", 0.0, ownedSurf(L, 0.0), "");
    w.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0) {
    auto s0 = Hdf5Reader::readSurfaceFilePerRank(stem, "water_depth", 0, size, NX * NY);
    REQUIRE(s0.extent(0) == static_cast<size_t>(NX * NY));
    for (int g = 0; g < NX * NY; ++g) {
      REQUIRE(s0(static_cast<size_t>(g)) == Approx(fsurf(g, 0.0)).margin(1e-9));
    }
  }
}

TEST_CASE("write a fixed file for the python/h5py interop check") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  IoLayout L = makeLayout(MPI_COMM_WORLD);
  const std::string path = std::string(kTmp) + "/h5py_check.h5";
  Grid grid(NX, NY, NZ, 80.0, 80.0, 0.1);
  Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::SerialGather);
  w.openFile(path, makeMeta());
  w.writeDomain(grid);
  w.writeSurfaceField("water_depth", 0.0, ownedSurf(L, 0.0), "");
  w.close();
  // No assertion here beyond not throwing; the h5py test (test_hdf5_h5py) validates content.
  REQUIRE(true);
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
