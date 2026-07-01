#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "frehg2_test.hpp"
#include "re/VanGenuchten.hpp"

using namespace frehg2;

TEST_CASE("VG: theta(h) and K(h) match closed forms") {
  SoilParams s;
  s.alpha = 1.43;
  s.n = 1.56;
  s.theta_s = 0.33;
  s.theta_r = 0.0;
  s.Ks_z = 2.89e-6;
  s.use_vg = true;
  s.use_mvg = false;
  s.air_entry = -0.02;

  const real h = -0.5;
  const real wc = VanGenuchten::waterContentFromHead(s, h);
  REQUIRE(wc > s.theta_r);
  REQUIRE(wc < s.theta_s);
  const real K = VanGenuchten::conductivityFromHead(s, h, s.Ks_z);
  REQUIRE(K > 0.0);
  REQUIRE(K <= s.Ks_z);
}

TEST_CASE("VG: wc -> h -> wc round-trip") {
  SoilParams s;
  s.alpha = 1.43;
  s.n = 1.56;
  s.theta_s = 0.33;
  s.theta_r = 0.0;
  s.use_vg = true;

  const std::vector<real> trials = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
  for (real wc0 : trials) {
    const real h = VanGenuchten::headFromWaterContent(s, wc0);
    const real wc1 = VanGenuchten::waterContentFromHead(s, h);
    REQUIRE(std::fabs(wc1 - wc0) < 1.0e-10);
  }
}

TEST_CASE("VG: MVG variant round-trip") {
  SoilParams s;
  s.alpha = 1.43;
  s.n = 1.56;
  s.theta_s = 0.33;
  s.theta_r = 0.0;
  s.use_vg = true;
  s.use_mvg = true;
  s.air_entry = -0.02;

  const real wc0 = 0.15;
  const real h = VanGenuchten::headFromWaterContent(s, wc0);
  const real wc1 = VanGenuchten::waterContentFromHead(s, h);
  REQUIRE(std::fabs(wc1 - wc0) < 1.0e-10);
}

TEST_CASE("VG: uniform wc=0.2 head consistency") {
  SoilParams s;
  s.alpha = 1.43;
  s.n = 1.56;
  s.theta_s = 0.33;
  s.theta_r = 0.0;
  s.use_vg = true;
  const real wc = 0.2;
  const real h = VanGenuchten::headFromWaterContent(s, wc);
  const real wc2 = VanGenuchten::waterContentFromHead(s, h);
  REQUIRE(std::fabs(wc2 - wc) < 1.0e-10);
}

int main() { return frehg2test::runAll(); }
