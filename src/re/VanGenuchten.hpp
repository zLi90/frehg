// Van Genuchten / modified VG constitutive relations (P5.1) — legacy-exact port of
// legacy/frehg/src/subroutines.c compute_wch/compute_hwc/compute_ch/compute_K/compute_dKdwc.
#ifndef FREHG2_RE_VAN_GENUCHTEN_HPP
#define FREHG2_RE_VAN_GENUCHTEN_HPP

#include <cmath>

#include <Kokkos_Core.hpp>

#include "frehg2/core/define.hpp"
#include "frehg2/re/SoilParams.hpp"

namespace frehg2 {

// All constitutive relations are KOKKOS_INLINE_FUNCTION so the P9 per-cell kernels that call
// them (RE flux/storage/corrector, coupling exchange) remain callable from device kernels in
// the P10 GPU pass. On the host backend these are ordinary inline calls -> bit-identical to P5.
struct VanGenuchten {
  KOKKOS_INLINE_FUNCTION static real m(real n) { return 1.0 - 1.0 / n; }

  KOKKOS_INLINE_FUNCTION static real wcm(const SoilParams& s) {
    const real mm = m(s.n);
    if (s.use_mvg) {
      return s.theta_r + (s.theta_s - s.theta_r) *
                         Kokkos::pow(1.0 + Kokkos::pow(Kokkos::fabs(s.air_entry) * s.alpha, s.n), mm);
    }
    return s.theta_s;
  }

  // theta(h) — legacy compute_wch
  KOKKOS_INLINE_FUNCTION static real waterContentFromHead(const SoilParams& s, real h) {
    real wc;
    if (s.use_vg) {
      const real mm = m(s.n);
      const real wcm_val = wcm(s);
      const real se = Kokkos::pow(1.0 + Kokkos::pow(Kokkos::fabs(s.alpha * h), s.n), -mm);
      if (h > s.air_entry)
        wc = s.theta_s;
      else
        wc = s.theta_r + (wcm_val - s.theta_r) * se;
    } else {
      wc = s.theta_r + (s.theta_s - s.theta_r) * Kokkos::exp(0.1634 * h);
    }
    if (wc > s.theta_s) wc = s.theta_s;
    if (wc < s.theta_r) wc = s.theta_r;
    return wc;
  }

  // h(theta) — legacy compute_hwc (head from water content stored in wc)
  KOKKOS_INLINE_FUNCTION static real headFromWaterContent(const SoilParams& s, real wc) {
    constexpr real eps = 1.0e-7;
    if (s.use_vg) {
      const real mm = m(s.n);
      const real wcm_val = wcm(s);
      real theta = wc;
      if (theta - s.theta_r < eps) theta = s.theta_r + eps;
      if (theta < s.theta_s) {
        return -(1.0 / s.alpha) *
               Kokkos::pow(Kokkos::pow((wcm_val - s.theta_r) / (theta - s.theta_r), 1.0 / mm) - 1.0,
                        1.0 / s.n);
      }
      return 0.0;
    }
    return Kokkos::log((wc - s.theta_r) / (s.theta_s - s.theta_r)) / 0.1634;
  }

  // d theta / d h — legacy compute_ch
  KOKKOS_INLINE_FUNCTION static real specificMoistureCapacity(const SoilParams& s, real h) {
    const real mm = m(s.n);
    const real wcm_val = wcm(s);
    const real nume = s.alpha * s.n * mm * (wcm_val - s.theta_r) *
                      Kokkos::pow(Kokkos::fabs(s.alpha * h), s.n - 1.0);
    const real deno = Kokkos::pow(1.0 + Kokkos::pow(Kokkos::fabs(s.alpha * h), s.n), mm + 1.0);
    real c = nume / deno;
    if (s.use_mvg) {
      if (h > s.air_entry) c = 0.0;
    } else if (h > 0.0) {
      c = 0.0;
    }
    return c;
  }

  // K(h) — legacy compute_K
  KOKKOS_INLINE_FUNCTION static real conductivityFromHead(const SoilParams& s, real h, real Ks) {
    if (!s.use_vg) return Ks * Kokkos::exp(0.1634 * h);
    const real mm = m(s.n);
    const real se = Kokkos::pow(1.0 + Kokkos::pow(Kokkos::fabs(s.alpha * h), s.n), -mm);
    const real wcm_val = wcm(s);
    real Keff;
    if (s.use_mvg) {
      const real nume =
          1.0 - Kokkos::pow(1.0 - Kokkos::pow(se * (s.theta_s - s.theta_r) / (wcm_val - s.theta_r),
                                         1.0 / mm),
                         mm);
      const real deno =
          1.0 - Kokkos::pow(1.0 - Kokkos::pow((s.theta_s - s.theta_r) / (wcm_val - s.theta_r),
                                          1.0 / mm),
                         mm);
      Keff = (deno == 0.0) ? Ks : Ks * Kokkos::pow(se, 0.5) * Kokkos::pow(nume / deno, 2.0);
    } else {
      Keff = Ks * Kokkos::pow(se, 0.5) *
             Kokkos::pow(1.0 - Kokkos::pow(1.0 - Kokkos::pow(se, 1.0 / mm), mm), 2.0);
    }
    if (Keff > Ks) Keff = Ks;
    if (s.use_mvg) {
      if (h > s.air_entry) Keff = Ks;
    } else if (h > 0.0) {
      Keff = Ks;
    }
    return Keff;
  }

  // dK/d theta — legacy compute_dKdwc (VG path, use_mvg==0 branch used by b2-gw)
  KOKKOS_INLINE_FUNCTION static real dKdwc(const SoilParams& s, real wc, real Ks) {
    real theta = wc;
    if (theta > 0.9999 * s.theta_s && theta < s.theta_s) theta = 0.9999 * s.theta_s;
    const real mm = m(s.n);
    const real se = (theta - s.theta_r) / (s.theta_s - s.theta_r);
    if (s.use_mvg) {
      const real wcm_val = wcm(s);
      const real c2 = (s.theta_s - s.theta_r) / (wcm_val - s.theta_r);
      const real c1 = 1.0 / Kokkos::pow(1.0 - Kokkos::pow(1.0 - Kokkos::pow(c2, 1.0 / mm), mm), 2.0);
      const real term0 = Kokkos::pow(1.0 - Kokkos::pow(c2 * se, 1.0 / mm), mm);
      const real term1 = 0.5 * c1 * Ks * Kokkos::pow(se, -0.5) * (1.0 - term0) * (1.0 - term0);
      const real term2 =
          2.0 * c1 * Ks * c2 * Kokkos::pow(se, 0.5) * Kokkos::pow(c2 * se, 1.0 / mm - 1.0) *
          (1.0 - term0) * Kokkos::pow(1.0 - Kokkos::pow(c2 * se, 1.0 / mm), mm - 1.0);
      return (term1 + term2) / (s.theta_s - s.theta_r);
    }
    const real term0 = Kokkos::pow(1.0 - Kokkos::pow(se, 1.0 / mm), mm);
    const real term1 = 0.5 * Ks * Kokkos::pow(se, -0.5) * (1.0 - term0) * (1.0 - term0);
    const real term2 =
        2.0 * Ks * Kokkos::pow(se, (2.0 - mm) / (2.0 * mm)) * (1.0 - term0) *
        Kokkos::pow(1.0 - Kokkos::pow(se, 1.0 / mm), mm - 1.0);
    return (term1 + term2) / (s.theta_s - s.theta_r);
  }
};

}  // namespace frehg2

#endif  // FREHG2_RE_VAN_GENUCHTEN_HPP
