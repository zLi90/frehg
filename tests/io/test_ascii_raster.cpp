// P3.3 acceptance: ESRI ASCII raster header + data + NODATA handling.
#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"
#include "io/AsciiRaster.hpp"

using namespace frehg2;

namespace {
// 3 cols x 2 rows; one NODATA cell. Lower-left CORNER variant.
const char* kCorner =
    "ncols 3\n"
    "nrows 2\n"
    "xllcorner 100.0\n"
    "yllcorner 200.0\n"
    "cellsize 10.0\n"
    "NODATA_value -9999\n"
    "1.0 2.0 3.0\n"
    "4.0 -9999 6.0\n";

// Same grid, lower-left CENTER variant (xll/yll shift by half a cell).
const char* kCenter =
    "ncols 3\n"
    "nrows 2\n"
    "xllcenter 105.0\n"
    "yllcenter 205.0\n"
    "cellsize 10.0\n"
    "NODATA_value -9999\n"
    "1.0 2.0 3.0\n"
    "4.0 -9999 6.0\n";
}  // namespace

TEST_CASE("header parsed (corner variant)") {
  AsciiRaster r = AsciiRaster::fromString(kCorner);
  REQUIRE(r.ncols() == 3);
  REQUIRE(r.nrows() == 2);
  REQUIRE(r.cellsize() == Approx(10.0).margin(1e-12));
  REQUIRE(r.xll() == Approx(100.0).margin(1e-9));
  REQUIRE(r.yll() == Approx(200.0).margin(1e-9));
  REQUIRE(r.noData() == Approx(-9999.0).margin(1e-9));
}

TEST_CASE("data values row-major (i inner, j outer)") {
  AsciiRaster r = AsciiRaster::fromString(kCorner);
  const auto& d = r.data();
  REQUIRE(d.extent(0) == 6u);
  REQUIRE(d(0) == Approx(1.0).margin(1e-12));  // (i=0,j=0)
  REQUIRE(d(2) == Approx(3.0).margin(1e-12));  // (i=2,j=0)
  REQUIRE(d(3) == Approx(4.0).margin(1e-12));  // (i=0,j=1)
  REQUIRE(d(5) == Approx(6.0).margin(1e-12));  // (i=2,j=1)
}

TEST_CASE("NODATA cell identified, valid cells not") {
  AsciiRaster r = AsciiRaster::fromString(kCorner);
  REQUIRE(r.isNoData(4, r.data()(4)));
  REQUIRE_FALSE(r.isNoData(0, r.data()(0)));
  auto mask = r.activeMask();
  REQUIRE(mask(0) == Approx(1.0).margin(1e-12));
  REQUIRE(mask(4) == Approx(0.0).margin(1e-12));
  REQUIRE(mask(5) == Approx(1.0).margin(1e-12));
}

TEST_CASE("xllcenter converted to corner consistently") {
  AsciiRaster r = AsciiRaster::fromString(kCenter);
  REQUIRE(r.xll() == Approx(100.0).margin(1e-9));  // 105 - 0.5*10
  REQUIRE(r.yll() == Approx(200.0).margin(1e-9));
  // data identical to corner variant
  REQUIRE(r.data()(0) == Approx(1.0).margin(1e-12));
  REQUIRE(r.data()(5) == Approx(6.0).margin(1e-12));
}
