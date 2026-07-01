#include "solute/Advection.hpp"

#include <algorithm>
#include <cmath>

#include "frehg2/core/ParallelFor.hpp"

namespace frehg2 {

namespace {

bool useMuscl(const SoluteParams& p) { return p.advection_scheme == "muscl"; }

KOKKOS_INLINE_FUNCTION real minmod(real a, real b) {
  if (a * b <= 0.0) return 0.0;
  return (Kokkos::fabs(a) < Kokkos::fabs(b)) ? a : b;
}

// Upwinded (optionally MUSCL-reconstructed) concentration at a face with velocity uf between
// left cell cL and right cell cR; cLL/cRR are the next cells out for the limiter slope.
KOKKOS_INLINE_FUNCTION real faceValue(bool muscl, real uf, real cLL, real cL, real cR,
                                      real cRR) {
  if (uf > 0.0) {
    if (!muscl) return cL;
    return cL + 0.5 * minmod(cL - cLL, cR - cL);
  }
  if (uf < 0.0) {
    if (!muscl) return cR;
    return cR - 0.5 * minmod(cRR - cR, cR - cL);
  }
  return 0.0;  // uf == 0: face carries no flux, value unused
}

}  // namespace

real advectSurface(RealArr1DHost& conc, const SoluteFlow& flow, const Grid& grid,
                   const SoluteParams& p, real dt) {
  if (!flow.hasSurface()) return 0.0;
  const int nx = grid.nx();
  const int ny = grid.ny();
  const real dx = grid.dx();
  const real dy = grid.dy();
  const Grid g = grid;
  RealArr1DHost u = flow.u;
  RealArr1DHost v = flow.v;

  // Max Courant number: sum of OUTGOING face fluxes per cell (the stability/positivity
  // limit of the explicit upwind update); interior faces only (outer boundary is zero-flux).
  // parallel_reduce with Kokkos::Max is order-independent -> bit-identical to the serial scan.
  real cfl = 0.0;
  Kokkos::parallel_reduce(
      "solute_adv_surf_cfl",
      Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<2>>({0, 0}, {nx, ny}),
      KOKKOS_LAMBDA(int i, int j, real& lmax) {
        const real uE = (i + 1 < nx) ? u(static_cast<size_t>(g.getSurfaceIndex(i, j))) : 0.0;
        const real uW = (i - 1 >= 0) ? u(static_cast<size_t>(g.getSurfaceIndex(i - 1, j))) : 0.0;
        const real vN = (j + 1 < ny) ? v(static_cast<size_t>(g.getSurfaceIndex(i, j))) : 0.0;
        const real vS = (j - 1 >= 0) ? v(static_cast<size_t>(g.getSurfaceIndex(i, j - 1))) : 0.0;
        const real out_x = Kokkos::max(uE, 0.0) + Kokkos::max(-uW, 0.0);
        const real out_y = Kokkos::max(vN, 0.0) + Kokkos::max(-vS, 0.0);
        const real c = out_x * dt / dx + out_y * dt / dy;
        if (c > lmax) lmax = c;
      },
      Kokkos::Max<real>(cfl));
  if (cfl > p.cfl_max) return cfl;  // refuse: leave concentration unchanged

  const bool muscl = useMuscl(p);
  RealArr1DHost cold("solute_cold_surf", static_cast<size_t>(nx) * static_cast<size_t>(ny));
  RealArr1DHost cc = conc;
  parallelForSurface("solute_adv_surf_copy", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    cold(static_cast<size_t>(i + j * nx)) = cc(static_cast<size_t>(g.getSurfaceIndex(i, j)));
  });

  parallelForSurface("solute_adv_surf_update", nx, ny, KOKKOS_LAMBDA(int i, int j) {
    auto C = [&](int ii, int jj) { return cold(static_cast<size_t>(ii + jj * nx)); };
    real Fe = 0.0, Fw = 0.0, Fn = 0.0, Fs = 0.0;
    if (i + 1 < nx) {  // east face: L=(i), R=(i+1)
      const real uf = u(static_cast<size_t>(g.getSurfaceIndex(i, j)));
      const real cLL = (i - 1 >= 0) ? C(i - 1, j) : C(i, j);
      const real cRR = (i + 2 < nx) ? C(i + 2, j) : C(i + 1, j);
      Fe = uf * faceValue(muscl, uf, cLL, C(i, j), C(i + 1, j), cRR);
    }
    if (i - 1 >= 0) {  // west face: L=(i-1), R=(i)
      const real uf = u(static_cast<size_t>(g.getSurfaceIndex(i - 1, j)));
      const real cLL = (i - 2 >= 0) ? C(i - 2, j) : C(i - 1, j);
      const real cRR = (i + 1 < nx) ? C(i + 1, j) : C(i, j);
      Fw = uf * faceValue(muscl, uf, cLL, C(i - 1, j), C(i, j), cRR);
    }
    if (j + 1 < ny) {  // north face
      const real uf = v(static_cast<size_t>(g.getSurfaceIndex(i, j)));
      const real cLL = (j - 1 >= 0) ? C(i, j - 1) : C(i, j);
      const real cRR = (j + 2 < ny) ? C(i, j + 2) : C(i, j + 1);
      Fn = uf * faceValue(muscl, uf, cLL, C(i, j), C(i, j + 1), cRR);
    }
    if (j - 1 >= 0) {  // south face
      const real uf = v(static_cast<size_t>(g.getSurfaceIndex(i, j - 1)));
      const real cLL = (j - 2 >= 0) ? C(i, j - 2) : C(i, j - 1);
      const real cRR = (j + 1 < ny) ? C(i, j + 1) : C(i, j);
      Fs = uf * faceValue(muscl, uf, cLL, C(i, j - 1), C(i, j), cRR);
    }
    cc(static_cast<size_t>(g.getSurfaceIndex(i, j))) =
        C(i, j) - dt / dx * (Fe - Fw) - dt / dy * (Fn - Fs);
  });
  return cfl;
}

real advectSubsurface(RealArr1DHost& conc, const SoluteFlow& flow, const Grid& grid,
                      const SoluteParams& p, real dt) {
  if (!flow.hasSubsurface()) return 0.0;
  const int nx = grid.nx();
  const int ny = grid.ny();
  const int nz = grid.nz();
  const real dx = grid.dx();
  const real dy = grid.dy();
  const real dz = grid.dz();
  const Grid g = grid;
  RealArr1DHost qx = flow.qx;
  RealArr1DHost qy = flow.qy;
  RealArr1DHost qz = flow.qz;

  real cfl = 0.0;
  Kokkos::parallel_reduce(
      "solute_adv_subs_cfl",
      Kokkos::MDRangePolicy<LoopExec, Kokkos::Rank<3>>({0, 0, 0}, {nx, ny, nz}),
      KOKKOS_LAMBDA(int i, int j, int k, real& lmax) {
        const real qE = (i + 1 < nx) ? qx(static_cast<size_t>(g.getIndex(i, j, k))) : 0.0;
        const real qW = (i - 1 >= 0) ? qx(static_cast<size_t>(g.getIndex(i - 1, j, k))) : 0.0;
        const real qN = (j + 1 < ny) ? qy(static_cast<size_t>(g.getIndex(i, j, k))) : 0.0;
        const real qS = (j - 1 >= 0) ? qy(static_cast<size_t>(g.getIndex(i, j - 1, k))) : 0.0;
        const real qD = (k + 1 < nz) ? qz(static_cast<size_t>(g.getIndex(i, j, k))) : 0.0;
        const real qU = (k - 1 >= 0) ? qz(static_cast<size_t>(g.getIndex(i, j, k - 1))) : 0.0;
        const real out_x = Kokkos::max(qE, 0.0) + Kokkos::max(-qW, 0.0);
        const real out_y = Kokkos::max(qN, 0.0) + Kokkos::max(-qS, 0.0);
        const real out_z = Kokkos::max(qD, 0.0) + Kokkos::max(-qU, 0.0);
        const real c = out_x * dt / dx + out_y * dt / dy + out_z * dt / dz;
        if (c > lmax) lmax = c;
      },
      Kokkos::Max<real>(cfl));
  if (cfl > p.cfl_max) return cfl;

  const bool muscl = useMuscl(p);
  const size_t n3 = static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
  RealArr1DHost cold("solute_cold_subs", n3);
  RealArr1DHost cc = conc;
  parallelForVolume("solute_adv_subs_copy", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    cold(static_cast<size_t>(i + j * nx + k * nx * ny)) =
        cc(static_cast<size_t>(g.getIndex(i, j, k)));
  });

  parallelForVolume("solute_adv_subs_update", nx, ny, nz, KOKKOS_LAMBDA(int i, int j, int k) {
    auto C = [&](int ii, int jj, int kk) {
      return cold(static_cast<size_t>(ii + jj * nx + kk * nx * ny));
    };
    real Fe = 0.0, Fw = 0.0, Fn = 0.0, Fs = 0.0, Fd = 0.0, Fu = 0.0;
    if (i + 1 < nx) {
      const real uf = qx(static_cast<size_t>(g.getIndex(i, j, k)));
      const real cLL = (i - 1 >= 0) ? C(i - 1, j, k) : C(i, j, k);
      const real cRR = (i + 2 < nx) ? C(i + 2, j, k) : C(i + 1, j, k);
      Fe = uf * faceValue(muscl, uf, cLL, C(i, j, k), C(i + 1, j, k), cRR);
    }
    if (i - 1 >= 0) {
      const real uf = qx(static_cast<size_t>(g.getIndex(i - 1, j, k)));
      const real cLL = (i - 2 >= 0) ? C(i - 2, j, k) : C(i - 1, j, k);
      const real cRR = (i + 1 < nx) ? C(i + 1, j, k) : C(i, j, k);
      Fw = uf * faceValue(muscl, uf, cLL, C(i - 1, j, k), C(i, j, k), cRR);
    }
    if (j + 1 < ny) {
      const real uf = qy(static_cast<size_t>(g.getIndex(i, j, k)));
      const real cLL = (j - 1 >= 0) ? C(i, j - 1, k) : C(i, j, k);
      const real cRR = (j + 2 < ny) ? C(i, j + 2, k) : C(i, j + 1, k);
      Fn = uf * faceValue(muscl, uf, cLL, C(i, j, k), C(i, j + 1, k), cRR);
    }
    if (j - 1 >= 0) {
      const real uf = qy(static_cast<size_t>(g.getIndex(i, j - 1, k)));
      const real cLL = (j - 2 >= 0) ? C(i, j - 2, k) : C(i, j - 1, k);
      const real cRR = (j + 1 < ny) ? C(i, j + 1, k) : C(i, j, k);
      Fs = uf * faceValue(muscl, uf, cLL, C(i, j - 1, k), C(i, j, k), cRR);
    }
    if (k + 1 < nz) {  // down face: L=(k), R=(k+1)
      const real uf = qz(static_cast<size_t>(g.getIndex(i, j, k)));
      const real cLL = (k - 1 >= 0) ? C(i, j, k - 1) : C(i, j, k);
      const real cRR = (k + 2 < nz) ? C(i, j, k + 2) : C(i, j, k + 1);
      Fd = uf * faceValue(muscl, uf, cLL, C(i, j, k), C(i, j, k + 1), cRR);
    }
    if (k - 1 >= 0) {  // up face: L=(k-1), R=(k)
      const real uf = qz(static_cast<size_t>(g.getIndex(i, j, k - 1)));
      const real cLL = (k - 2 >= 0) ? C(i, j, k - 2) : C(i, j, k - 1);
      const real cRR = (k + 1 < nz) ? C(i, j, k + 1) : C(i, j, k);
      Fu = uf * faceValue(muscl, uf, cLL, C(i, j, k - 1), C(i, j, k), cRR);
    }
    cc(static_cast<size_t>(g.getIndex(i, j, k))) = C(i, j, k) - dt / dx * (Fe - Fw) -
                                                   dt / dy * (Fn - Fs) - dt / dz * (Fd - Fu);
  });
  return cfl;
}

}  // namespace frehg2
