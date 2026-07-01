// P13.3.2 acceptance: the soil YAML schema parses 2-class, 4-class, and 16-class `soil.types`
// lists into SoilParams tuples (per-class VG/Ksat; global groundwater.* defaults + per-class
// overrides). Also asserts fail-loud on a class missing a required field.
#define FREHG2_TEST_IMPL
#include <string>

#include "frehg2_test.hpp"
#include "io/Config.hpp"
#include "soil/SoilConfig.hpp"

using namespace frehg2;

namespace {
// Build a config string with `n` soil classes, each with a distinct Ksat and theta_s so we can
// verify ordered, per-class parsing.
std::string makeCfg(int n) {
  std::string s =
      "schema_version: '2.0'\n"
      "simulation: {id: soil, mode: groundwater}\n"
      "domain: {nx: 4, ny: 4, nz: 4, dx: 1.0, dy: 1.0, dz: 0.1}\n"
      "time: {dt: 1.0, t_end: 1.0}\n"
      "modules: {surface_water: false, groundwater: true, solute: false}\n"
      "groundwater: {specific_storage: 2.0e-5, use_vg: true, air_entry_value: -0.03}\n"
      "soil:\n  types:\n";
  for (int c = 0; c < n; ++c) {
    const double ks = 1.0e-6 * static_cast<double>(c + 1);
    const double ts = 0.30 + 0.01 * static_cast<double>(c);
    s += "    - {id: " + std::to_string(c) + ", theta_s: " + std::to_string(ts) +
         ", theta_r: 0.05, vg: {alpha: 1.4, n: 1.5}, k_sat: {x: " + std::to_string(ks) +
         ", y: " + std::to_string(ks) + ", z: " + std::to_string(ks) + "}}\n";
  }
  return s;
}
}  // namespace

TEST_CASE("yaml soil: parse 2 / 4 / 16 classes, ordered and complete") {
  for (int n : {2, 4, 16}) {
    Config cfg = Config::fromString(makeCfg(n), "");
    std::vector<SoilParams> classes = parseSoilClasses(cfg);
    REQUIRE(static_cast<int>(classes.size()) == n);
    for (int c = 0; c < n; ++c) {
      REQUIRE(classes[static_cast<size_t>(c)].Ks_z ==
              Approx(1.0e-6 * (c + 1)).margin(1e-15));
      REQUIRE(classes[static_cast<size_t>(c)].theta_s ==
              Approx(0.30 + 0.01 * c).margin(1e-6));
      REQUIRE(classes[static_cast<size_t>(c)].theta_r == Approx(0.05).margin(1e-12));
      REQUIRE(classes[static_cast<size_t>(c)].alpha == Approx(1.4).margin(1e-12));
      REQUIRE(classes[static_cast<size_t>(c)].n == Approx(1.5).margin(1e-12));
      // Global groundwater.* defaults applied to every class.
      REQUIRE(classes[static_cast<size_t>(c)].Ss == Approx(2.0e-5).margin(1e-12));
      REQUIRE(classes[static_cast<size_t>(c)].air_entry == Approx(-0.03).margin(1e-12));
      REQUIRE(classes[static_cast<size_t>(c)].use_vg);
    }
  }
}

TEST_CASE("yaml soil: per-class override of storage / use_vg") {
  const char* cfg = R"(schema_version: '2.0'
simulation: {id: soil, mode: groundwater}
domain: {nx: 2, ny: 2, nz: 2, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: false, groundwater: true, solute: false}
groundwater: {specific_storage: 1.0e-5}
soil:
  types:
    - {id: 0, theta_s: 0.40, theta_r: 0.0, vg: {alpha: 1.4, n: 1.5}, k_sat: {x: 1.0e-6, y: 1.0e-6, z: 1.0e-6}}
    - {id: 1, theta_s: 0.35, theta_r: 0.0, vg: {alpha: 2.0, n: 1.8}, k_sat: {x: 5.0e-6, y: 5.0e-6, z: 5.0e-6}, specific_storage: 9.0e-5, use_vg: false}
)";
  std::vector<SoilParams> classes = parseSoilClasses(Config::fromString(cfg, ""));
  REQUIRE(classes.size() == 2u);
  REQUIRE(classes[0].Ss == Approx(1.0e-5).margin(1e-12));  // global default
  REQUIRE(classes[0].use_vg);
  REQUIRE(classes[1].Ss == Approx(9.0e-5).margin(1e-12));  // per-class override
  REQUIRE_FALSE(classes[1].use_vg);
}

TEST_CASE("yaml soil: missing field and empty list fail loud") {
  const char* bad = R"(schema_version: '2.0'
simulation: {id: soil, mode: groundwater}
domain: {nx: 2, ny: 2, nz: 2, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: false, groundwater: true, solute: false}
soil:
  types:
    - {id: 0, theta_r: 0.0, vg: {alpha: 1.4, n: 1.5}, k_sat: {x: 1.0e-6, y: 1.0e-6, z: 1.0e-6}}
)";
  bool threw = false;
  try {
    parseSoilClasses(Config::fromString(bad, ""));
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);

  const char* empty = R"(schema_version: '2.0'
simulation: {id: soil, mode: groundwater}
domain: {nx: 2, ny: 2, nz: 2, dx: 1.0, dy: 1.0, dz: 0.1}
time: {dt: 1.0, t_end: 1.0}
modules: {surface_water: false, groundwater: true, solute: false}
)";
  threw = false;
  try {
    parseSoilClasses(Config::fromString(empty, ""));
  } catch (const std::exception&) {
    threw = true;
  }
  REQUIRE(threw);
}
