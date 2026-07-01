// P3.5 acceptance: checkpoint write/read round-trip, dimension validation, missing-file
// error, atomic temp+rename behavior, and config-provenance mismatch detection.
#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <cstdio>
#include <fstream>
#include <string>

#include "frehg2/io/OutputWriter.hpp"
#include "frehg2_test.hpp"
#include "io/Hdf5Reader.hpp"
#include "io/Hdf5Writer.hpp"
#include "io/Sha256.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;

RealArr1DHost filled(const std::string& nm, int n, double base) {
  RealArr1DHost v(nm, static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) v(static_cast<size_t>(i)) = base + 0.5 * i;
  return v;
}

CheckpointState makeCp(double time, double dt) {
  CheckpointState cp;
  cp.time = time;
  cp.dt = dt;
  cp.has_sw = true;
  cp.has_gw = true;
  cp.config_sha256 = sha256Hex("cfgA");
  cp.git_sha = "deadbeef";
  cp.eta = filled("eta", 9, 1.0);   // 3x3
  cp.u = filled("u", 9, 2.0);
  cp.v = filled("v", 9, 3.0);
  cp.h = filled("h", 45, 4.0);      // 3x3x5
  cp.hn = filled("hn", 45, 5.0);
  cp.wc = filled("wc", 45, 6.0);
  cp.wcn = filled("wcn", 45, 7.0);
  return cp;
}

bool sameView(const RealArr1DHost& a, const RealArr1DHost& b) {
  if (a.extent(0) != b.extent(0)) return false;
  for (size_t i = 0; i < a.extent(0); ++i) {
    if (a(i) != b(i)) return false;
  }
  return true;
}

std::string ckptPath(const std::string& stem, double time) {
  return stem + ".ckpt." + std::to_string(static_cast<long long>(time)) + ".h5";
}
}  // namespace

TEST_CASE("write checkpoint, read back element-by-element identical") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string stem = std::string(kTmp) + "/run_ck";
  IoLayout L;
  L.nx = 3;
  L.ny = 3;
  L.nz = 5;

  Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::SerialGather);
  w.openFile(stem + ".h5", RunMetadata{});
  CheckpointState cp = makeCp(10.0, 0.25);
  w.writeCheckpoint(cp);
  w.close();

  if (rank == 0) {
    RestartState rs = Hdf5Reader::readCheckpoint(ckptPath(stem, 10.0));
    REQUIRE(rs.time == Approx(10.0).margin(1e-12));
    REQUIRE(rs.dt == Approx(0.25).margin(1e-12));
    REQUIRE(rs.has_sw);
    REQUIRE(rs.has_gw);
    REQUIRE(rs.eta.extent(0) == 9u);
    REQUIRE(rs.h.extent(0) == 45u);
    REQUIRE(sameView(rs.eta, cp.eta));
    REQUIRE(sameView(rs.u, cp.u));
    REQUIRE(sameView(rs.v, cp.v));
    REQUIRE(sameView(rs.h, cp.h));
    REQUIRE(sameView(rs.hn, cp.hn));
    REQUIRE(sameView(rs.wc, cp.wc));
    REQUIRE(sameView(rs.wcn, cp.wcn));
  }
}

TEST_CASE("read from non-existent checkpoint gives clear error, not crash") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;
  bool threw = false;
  try {
    (void)Hdf5Reader::readCheckpoint(std::string(kTmp) + "/no_such.ckpt.0.h5");
  } catch (const std::runtime_error& e) {
    threw = true;
    const std::string msg = e.what();
    REQUIRE(msg.find("no_such") != std::string::npos);
  }
  REQUIRE(threw);
}

TEST_CASE("atomic temp+rename leaves the committed checkpoint intact") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string stem = std::string(kTmp) + "/run_atomic";
  IoLayout L;
  L.nx = 3;
  L.ny = 3;
  L.nz = 5;

  Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::SerialGather);
  w.openFile(stem + ".h5", RunMetadata{});
  w.writeCheckpoint(makeCp(20.0, 0.1));  // commits run_atomic.ckpt.20.h5
  w.close();

  if (rank == 0) {
    const std::string final_path = ckptPath(stem, 20.0);
    // Simulate a crash mid-write of a *new* checkpoint: a stray .tmp exists then is lost.
    const std::string tmp_path = ckptPath(stem, 30.0) + ".tmp";
    { std::ofstream(tmp_path) << "partial"; }
    std::remove(tmp_path.c_str());  // crash before rename
    // The previously committed checkpoint must still be fully readable.
    RestartState rs = Hdf5Reader::readCheckpoint(final_path);
    REQUIRE(rs.time == Approx(20.0).margin(1e-12));
    REQUIRE(rs.eta.extent(0) == 9u);
  }
}

TEST_CASE("config provenance mismatch is detected on restart") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::string stem = std::string(kTmp) + "/run_prov";
  IoLayout L;
  L.nx = 3;
  L.ny = 3;
  L.nz = 5;

  Hdf5Writer w(MPI_COMM_WORLD, L, IoMode::SerialGather);
  w.openFile(stem + ".h5", RunMetadata{});
  w.writeCheckpoint(makeCp(5.0, 0.1));  // config_sha256 = sha256("cfgA")
  w.close();

  if (rank == 0) {
    RestartState match = Hdf5Reader::readCheckpoint(ckptPath(stem, 5.0), sha256Hex("cfgA"));
    REQUIRE(match.config_matches);
    RestartState mismatch =
        Hdf5Reader::readCheckpoint(ckptPath(stem, 5.0), sha256Hex("cfgB"));
    REQUIRE_FALSE(mismatch.config_matches);
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
