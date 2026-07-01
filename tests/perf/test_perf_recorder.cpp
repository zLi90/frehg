// P21 unit test: PerfRecorder accumulates named regions and counters with numerical bounds.
#include "frehg2/perf/Counters.hpp"
#include "frehg2/perf/PerfRecorder.hpp"
#include "frehg2/perf/Timer.hpp"
#include "frehg2_test.hpp"

#include <chrono>
#include <sstream>
#include <thread>

using namespace frehg2;
using namespace frehg2::perf;

TEST_CASE("PerfRecorder accumulates region times and counters") {
  PerfRecorder rec;
  rec.addTime(Region::SwAssembly, 0.12);
  rec.addTime(Region::SwKsp, 0.34);
  rec.counters().addCells(100);
  rec.counters().addKspIters(7);
  rec.counters().addBytes(4096);

  REQUIRE(rec.timeSeconds(Region::SwAssembly) == Approx(0.12).margin(1e-15));
  REQUIRE(rec.timeSeconds(Region::SwKsp) == Approx(0.34).margin(1e-15));
  REQUIRE(rec.counters().cells_touched == 100);
  REQUIRE(rec.counters().ksp_iterations == 7);
  REQUIRE(rec.counters().bytes_staged == 4096);

  std::ostringstream ss;
  rec.writeSummaryLines(ss, 10, 0);
  const std::string text = ss.str();
  REQUIRE(text.find("perf_sw_assembly_seconds 0.12") != std::string::npos);
  REQUIRE(text.find("perf_sw_ksp_seconds 0.34") != std::string::npos);
  REQUIRE(text.find("perf_cells_touched 100") != std::string::npos);
  REQUIRE(text.find("perf_ksp_iterations 7") != std::string::npos);
}

TEST_CASE("ScopedTimer records elapsed time on destruction") {
  PerfRecorder rec;
  {
    ScopedTimer t(&rec, Region::Io);
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
  }
  REQUIRE(rec.timeSeconds(Region::Io) > 0.001);
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  return frehg2test::runAll();
}
