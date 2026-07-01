// P13.3.1 acceptance: SoilMap data structure — class table + per-cell class index, uniform path,
// CSV load (1k-1M rows), classAt lookup, and fail-loud validation.
#define FREHG2_TEST_IMPL
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "frehg2_test.hpp"
#include "soil/SoilMap.hpp"

using namespace frehg2;

namespace {
SoilParams makeClass(double ks, double ts) {
  SoilParams sp;
  sp.Ks_z = ks;
  sp.Ks_x = ks;
  sp.Ks_y = ks;
  sp.theta_s = ts;
  sp.theta_r = 0.05;
  return sp;
}

std::string writeTemp(const std::string& name, const std::string& body) {
  std::string path = std::string("/tmp/") + name;
  std::ofstream o(path);
  o << body;
  o.close();
  return path;
}
}  // namespace

TEST_CASE("soil map: uniform single class reproduces every cell") {
  SoilMap m;
  m.setUniform(8, 5, makeClass(1.0e-5, 0.40));
  REQUIRE(m.nClasses() == 1);
  REQUIRE(m.nx() == 8);
  REQUIRE(m.ny() == 5);
  REQUIRE(m.isUniform());
  for (int j = 0; j < 5; ++j)
    for (int i = 0; i < 8; ++i) {
      REQUIRE(m.classIdAt(i, j) == 0);
      REQUIRE(m.classAt(i, j).Ks_z == Approx(1.0e-5).margin(1e-18));
      REQUIRE(m.classAt(i, j).theta_s == Approx(0.40).margin(1e-15));
    }
}

TEST_CASE("soil map: explicit multi-class index dispatches per cell") {
  SoilMap m;
  m.setClasses({makeClass(1.0e-5, 0.40), makeClass(5.0e-6, 0.30), makeClass(2.0e-5, 0.45)});
  // 3x2 grid (row-major i + j*nx): left column class 0, middle class 1, right class 2.
  std::vector<int> idx = {0, 1, 2,
                          0, 1, 2};
  m.setClassIndex(3, 2, idx);
  REQUIRE(m.nClasses() == 3);
  REQUIRE_FALSE(m.isUniform());
  REQUIRE(m.classIdAt(0, 0) == 0);
  REQUIRE(m.classIdAt(1, 1) == 1);
  REQUIRE(m.classIdAt(2, 0) == 2);
  REQUIRE(m.classAt(1, 0).theta_s == Approx(0.30).margin(1e-15));
  REQUIRE(m.classAt(2, 1).Ks_z == Approx(2.0e-5).margin(1e-18));
  REQUIRE(m.classById(0).Ks_z == Approx(1.0e-5).margin(1e-18));
}

TEST_CASE("soil map: fully heterogeneous 3D index dispatches per cell (P23)") {
  SoilMap m;
  m.setClasses({makeClass(1.0e-5, 0.40), makeClass(5.0e-6, 0.30)});
  // 2x2x3 grid, row-major i + j*nx + k*nx*ny: top layer (k=0) all class 0, then class 1, then
  // a per-cell checkerboard at k=2.
  const int nx = 2, ny = 2, nz = 3;
  std::vector<int> idx(static_cast<size_t>(nx * ny * nz), 0);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      idx[static_cast<size_t>(i + j * nx + 0 * nx * ny)] = 0;
      idx[static_cast<size_t>(i + j * nx + 1 * nx * ny)] = 1;
      idx[static_cast<size_t>(i + j * nx + 2 * nx * ny)] = (i + j) % 2;
    }
  m.setClassIndex3D(nx, ny, nz, idx);
  REQUIRE(m.is3D());
  REQUIRE(m.nz() == 3);
  REQUIRE(m.classIdAt(0, 0, 0) == 0);
  REQUIRE(m.classIdAt(1, 1, 1) == 1);
  REQUIRE(m.classIdAt(1, 0, 2) == 1);  // i+j odd
  REQUIRE(m.classIdAt(0, 0, 2) == 0);  // i+j even
  REQUIRE(m.classAt(1, 1, 1).theta_s == Approx(0.30).margin(1e-15));
  REQUIRE(m.classAt(0, 0, 0).theta_s == Approx(0.40).margin(1e-15));

  // A 2D (per-column) map ignores k => column class replicated across layers.
  SoilMap col;
  col.setClasses({makeClass(1.0e-5, 0.40), makeClass(5.0e-6, 0.30)});
  col.setClassIndex(2, 2, {0, 1, 1, 0});
  REQUIRE_FALSE(col.is3D());
  REQUIRE(col.classIdAt(1, 0, 0) == 1);
  REQUIRE(col.classIdAt(1, 0, 7) == 1);  // any k => same column class

  // 3D validation: wrong count and out-of-range id both throw.
  bool threw = false;
  try {
    m.setClassIndex3D(2, 2, 3, {0, 1});
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);
  threw = false;
  try {
    m.setClassIndex3D(1, 1, 1, {5});
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);
}

TEST_CASE("soil map: CSV load handles whitespace and comma, large grids") {
  SoilMap m;
  m.setClasses({makeClass(1.0e-5, 0.40), makeClass(5.0e-6, 0.30)});
  // 4x3 with mixed separators.
  const std::string body = "0,1,0,1\n1 0 1 0\n0,0,1,1\n";
  const std::string path = writeTemp("frehg2_soil_idx.csv", body);
  m.loadIndexFromCSV(path, 4, 3);
  REQUIRE(m.classIdAt(0, 0) == 0);
  REQUIRE(m.classIdAt(1, 0) == 1);
  REQUIRE(m.classIdAt(0, 1) == 1);
  REQUIRE(m.classIdAt(3, 2) == 1);
  std::remove(path.c_str());

  // Large grid: 1000x1000 = 1e6 cells, alternating 0/1.
  SoilMap big;
  big.setClasses({makeClass(1.0e-5, 0.40), makeClass(5.0e-6, 0.30)});
  const int N = 1000;
  std::string big_body;
  big_body.reserve(static_cast<size_t>(N) * N * 2);
  for (int k = 0; k < N * N; ++k) {
    big_body += (k % 2 == 0) ? '0' : '1';
    big_body += (k % N == N - 1) ? '\n' : ' ';
  }
  const std::string big_path = writeTemp("frehg2_soil_big.csv", big_body);
  big.loadIndexFromCSV(big_path, N, N);
  REQUIRE(big.nx() == N);
  REQUIRE(big.ny() == N);
  REQUIRE(big.classIdAt(0, 0) == 0);
  REQUIRE(big.classIdAt(1, 0) == 1);
  std::remove(big_path.c_str());
}

TEST_CASE("soil map: validation fails loud") {
  SoilMap m;
  // No classes set yet.
  bool threw = false;
  try {
    m.setClassIndex(2, 2, {0, 0, 0, 0});
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);

  m.setClasses({makeClass(1.0e-5, 0.40), makeClass(5.0e-6, 0.30)});

  // Wrong count.
  threw = false;
  try {
    m.setClassIndex(2, 2, {0, 1, 0});
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);

  // Out-of-range class id (only 2 classes => valid ids 0,1).
  threw = false;
  try {
    m.setClassIndex(2, 2, {0, 1, 2, 0});
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);

  // Empty class table.
  threw = false;
  try {
    SoilMap m2;
    m2.setClasses({});
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);
}
