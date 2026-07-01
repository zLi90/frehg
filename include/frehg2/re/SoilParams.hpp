// Soil hydraulic parameters for the Richards-equation constitutive relations (P5/P13).
//
// Promoted to a public header (from src/re/VanGenuchten.hpp) in P13 so the non-uniform
// `SoilMap` (src/soil/) can hold a vector of these per soil class WITHOUT creating a
// re <-> soil dependency cycle. The struct and its defaults are unchanged from P5.
#ifndef FREHG2_RE_SOIL_PARAMS_HPP
#define FREHG2_RE_SOIL_PARAMS_HPP

#include "frehg2/core/define.hpp"

namespace frehg2 {

struct SoilParams {
  real alpha = 1.43;     // VG alpha (legacy vga)
  real n = 1.56;         // VG n (legacy vgn)
  real theta_s = 0.33;   // wcs (saturated water content / porosity)
  real theta_r = 0.0;    // wcr (residual water content)
  real Ks_x = 0.0;
  real Ks_y = 0.0;
  real Ks_z = 2.89e-6;
  real Ss = 1.0e-5;        // specific storage
  bool use_vg = true;
  bool use_mvg = false;
  real air_entry = -0.02;  // legacy aev
};

}  // namespace frehg2

#endif  // FREHG2_RE_SOIL_PARAMS_HPP
