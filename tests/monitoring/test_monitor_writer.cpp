#define FREHG2_TEST_IMPL
#define FREHG2_TEST_USE_KOKKOS
#include "frehg2_test.hpp"

#include <fstream>
#include <sstream>
#include <string>

#include "core/Grid.hpp"
#include "monitoring/MonitorSpec.hpp"
#include "monitoring/MonitorWriter.hpp"
#include "swe/SweSolver.hpp"

using namespace frehg2;

namespace {
const char* kTmp = FREHG2_IO_TMP;

std::vector<std::string> splitCsvLine(const std::string& line) {
  std::vector<std::string> out;
  std::string cell;
  std::istringstream ss(line);
  while (std::getline(ss, cell, ',')) out.push_back(cell);
  return out;
}
}  // namespace

TEST_CASE("MonitorWriter CSV schema and probe values") {
  Grid g(5, 1, 1, 1.0, 1.0, 0.1);
  SweSolver swe(g, nullptr);
  swe.setParams(SweParams{});
  RealArr1DHost bed("bed", 5);
  for (size_t i = 0; i < 5; ++i) bed(i) = 0.0;
  swe.setBathymetry(bed);
  RealArr1DHost eta("eta", 5);
  for (size_t i = 0; i < 5; ++i) eta(i) = 2.0 + 0.1 * static_cast<double>(i);
  swe.setInitialEtaPhysical(eta);
  swe.setInitialVelocityConstant(1.0, 0.0);
  swe.finalizeInitialState();

  MonitorBundle bundle;
  ProbeSpec pr;
  pr.name = "mid";
  pr.indices_from_grid = true;
  pr.gi = 2;
  pr.gj = 0;
  pr.fields = {"eta", "depth", "u"};
  bundle.probes.push_back(pr);

  LineFluxSpec line;
  line.name = "outlet";
  line.p0 = {2.0, 0.5, 0.0};
  line.p1 = {4.0, 0.5, 0.0};
  line.field = "u";
  bundle.lines.push_back(line);

  const std::string csv_path = std::string(kTmp) + "/monitors/monitor_unit.csv";
  MonitorWriter writer;
  writer.configure(bundle, kTmp, "monitor_unit", 5, 1, 1, 1.0, 1.0, 0.1, 0.0, 0.0, 0.0,
                   nullptr);
  writer.open(false);
  writer.writeRow(0.0, &swe, nullptr, g, g);
  writer.writeRow(1.0, &swe, nullptr, g, g);
  writer.close();

  std::ifstream in(csv_path);
  REQUIRE(in.good());
  std::string header;
  REQUIRE(std::getline(in, header));
  REQUIRE(header ==
            "time,mid.eta,mid.depth,mid.u,outlet.flux");

  std::string row0;
  REQUIRE(std::getline(in, row0));
  const auto cells = splitCsvLine(row0);
  REQUIRE(cells.size() == 5);
  REQUIRE(std::stod(cells[0]) == Approx(0.0).margin(1e-12));
  REQUIRE(std::stod(cells[1]) == Approx(2.2).margin(1e-12));
  REQUIRE(std::stod(cells[2]) == Approx(2.2).margin(1e-12));
  REQUIRE(std::stod(cells[3]) == Approx(1.0).margin(1e-12));
  REQUIRE(std::stod(cells[4]) == Approx(2.0).margin(1e-6));
}
