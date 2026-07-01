// P8.3.6: the YAML `solute:` block parses into SoluteParams with the documented defaults,
// and explicit values override them. A minimal config omitting solute keys yields defaults.
#define FREHG2_TEST_IMPL
#include "frehg2_test.hpp"

#include <string>

#include "io/Config.hpp"
#include "solute/SoluteParams.hpp"

using namespace frehg2;

namespace {

const char* kMinimal = R"(
modules:
  shallow_water: true
time:
  t_end: 1.0
  dt: 0.1
domain:
  nx: 2
  ny: 2
  nz: 1
solute:
  enabled: true
)";

const char* kExplicit = R"(
modules:
  shallow_water: true
time:
  t_end: 1.0
  dt: 0.1
domain:
  nx: 2
  ny: 2
  nz: 1
solute:
  enabled: true
  c_rain: 12.5
  k_decay: 0.02
  D: 1.0e-6
  advection_scheme: muscl
  diffusion_scheme: implicit
  cfl_max: 0.8
)";

}  // namespace

TEST_CASE("yaml solute: defaults applied when keys omitted") {
  Config cfg = Config::fromString(kMinimal);
  SoluteParams p = SoluteParams::fromConfig(cfg);
  REQUIRE(p.enabled == true);
  REQUIRE(p.c_rain == Approx(0.0).margin(0.0));
  REQUIRE(p.k_decay == Approx(0.0).margin(0.0));
  REQUIRE(p.D == Approx(1.0e-9).epsilon(1e-12));
  REQUIRE(p.advection_scheme == std::string("upwind"));
  REQUIRE(p.diffusion_scheme == std::string("implicit"));
  REQUIRE(p.cfl_max == Approx(1.0).margin(0.0));
}

TEST_CASE("yaml solute: explicit values override defaults") {
  Config cfg = Config::fromString(kExplicit);
  SoluteParams p = SoluteParams::fromConfig(cfg);
  REQUIRE(p.enabled == true);
  REQUIRE(p.c_rain == Approx(12.5).margin(1e-12));
  REQUIRE(p.k_decay == Approx(0.02).margin(1e-12));
  REQUIRE(p.D == Approx(1.0e-6).epsilon(1e-9));
  REQUIRE(p.advection_scheme == std::string("muscl"));
  REQUIRE(p.cfl_max == Approx(0.8).margin(1e-12));
}

TEST_CASE("yaml solute: disabled by default when block absent") {
  const char* no_solute = R"(
modules:
  shallow_water: true
time:
  t_end: 1.0
  dt: 0.1
domain:
  nx: 2
  ny: 2
  nz: 1
)";
  Config cfg = Config::fromString(no_solute);
  SoluteParams p = SoluteParams::fromConfig(cfg);
  REQUIRE(p.enabled == false);
}
