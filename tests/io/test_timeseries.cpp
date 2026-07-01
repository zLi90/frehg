// P3.4 acceptance: time-series CSV/whitespace parsing + interpolation + clamping.
#define FREHG2_TEST_IMPL
#include "frehg2_test.hpp"
#include "io/TimeSeries.hpp"

using namespace frehg2;

TEST_CASE("CSV: interpolation at data points and midpoint") {
  TimeSeries ts = TimeSeries::fromString("0.0, 10.0\n1.0, 20.0\n2.0, 40.0\n");
  REQUIRE(ts.size() == 3u);
  REQUIRE(ts.getValueAt(0.0) == Approx(10.0).margin(1e-12));
  REQUIRE(ts.getValueAt(1.0) == Approx(20.0).margin(1e-12));
  REQUIRE(ts.getValueAt(2.0) == Approx(40.0).margin(1e-12));
  REQUIRE(ts.getValueAt(0.5) == Approx(15.0).margin(1e-12));
  REQUIRE(ts.getValueAt(1.5) == Approx(30.0).margin(1e-12));
}

TEST_CASE("space-separated and comments parse equivalently") {
  TimeSeries ts = TimeSeries::fromString("# t v\n0 10\n1 20\n\n2 40\n");
  REQUIRE(ts.size() == 3u);
  REQUIRE(ts.getValueAt(0.5) == Approx(15.0).margin(1e-12));
}

TEST_CASE("out-of-bounds clamps to nearest boundary value") {
  TimeSeries ts = TimeSeries::fromString("0 10\n1 20\n2 40\n");
  REQUIRE(ts.getValueAt(-5.0) == Approx(10.0).margin(1e-12));
  REQUIRE(ts.getValueAt(99.0) == Approx(40.0).margin(1e-12));
}

TEST_CASE("unsorted input is sorted before interpolation") {
  TimeSeries ts = TimeSeries::fromString("2 40\n0 10\n1 20\n");
  REQUIRE(ts.times().front() == Approx(0.0).margin(1e-12));
  REQUIRE(ts.times().back() == Approx(2.0).margin(1e-12));
  REQUIRE(ts.getValueAt(1.5) == Approx(30.0).margin(1e-12));
}

TEST_CASE("empty input throws") {
  bool threw = false;
  try {
    TimeSeries::fromString("# only a comment\n");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  REQUIRE(threw);
}
