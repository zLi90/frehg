// Halo-padded surface-water work arrays for the semi-implicit SWE solver (P4).
//
// The SWE solver owns its state in a flat, halo-padded layout (one ghost cell on each
// horizontal side), matching legacy Frehg's per-cell computations but using Frehg2's
// `Grid::getSurfaceIndex(i,j)` ordering rather than the legacy interior-first map. All
// arrays are length `Grid::nSurfaceStorageCell() == (nx+2)*(ny+2)`. Storage is host-side:
// SWE local loops stay plain serial/MPI through Phase 4 (DP3 defers on-node Kokkos/GPU
// parallelization to the global P9/P10 passes). Names mirror the legacy `Data` struct.
#ifndef FREHG2_SWE_SWE_FIELDS_HPP
#define FREHG2_SWE_SWE_FIELDS_HPP

#include <string>
#include <utility>
#include <vector>

#include "core/Grid.hpp"
#include "frehg2/core/define.hpp"

namespace frehg2 {

struct SweFields {
  explicit SweFields(const Grid& grid) {
    const size_t n = static_cast<size_t>(grid.nSurfaceStorageCell());
    auto a = [n](const char* name) { return RealArr1DHost(name, n); };
    eta = a("eta");
    etan = a("etan");
    bottom = a("bottom");
    bottomXP = a("bottomXP");
    bottomYP = a("bottomYP");
    dept = a("dept");
    deptx = a("deptx");
    depty = a("depty");
    uu = a("uu");
    vv = a("vv");
    un = a("un");
    vn = a("vn");
    uy = a("uy");
    vx = a("vx");
    Vs = a("Vs");
    Vsn = a("Vsn");
    Asx = a("Asx");
    Asy = a("Asy");
    Asz = a("Asz");
    Aszx = a("Aszx");
    Aszy = a("Aszy");
    Vsx = a("Vsx");
    Vsy = a("Vsy");
    CDx = a("CDx");
    CDy = a("CDy");
    Dx = a("Dx");
    Dy = a("Dy");
    Ex = a("Ex");
    Ey = a("Ey");
    cflx = a("cflx");
    cfly = a("cfly");
    Fu = a("Fu");
    Fv = a("Fv");
    rain = a("rain");
    evap = a("evap");
    Sxp = a("Sxp");
    Sxm = a("Sxm");
    Syp = a("Syp");
    Sym = a("Sym");
    Sct = a("Sct");
    Srhs = a("Srhs");
    cfl_active = a("cfl_active");
    wtfx = a("wtfx");
    wtfy = a("wtfy");
  }

  // Free surface and bathymetry.
  RealArr1DHost eta;       // water surface elevation (primary unknown)
  RealArr1DHost etan;      // eta at previous step
  RealArr1DHost bottom;    // bed elevation (cell)
  RealArr1DHost bottomXP;  // bed elevation at x+ face (max of neighbors)
  RealArr1DHost bottomYP;  // bed elevation at y+ face

  // Depths.
  RealArr1DHost dept;   // cell depth = eta - bottom (>= 0)
  RealArr1DHost deptx;  // x+ face depth
  RealArr1DHost depty;  // y+ face depth

  // Velocities (face) and corner-interpolated velocities.
  RealArr1DHost uu;  // x face velocity
  RealArr1DHost vv;  // y face velocity
  RealArr1DHost un;  // previous uu
  RealArr1DHost vn;  // previous vv
  RealArr1DHost uy;  // u interpolated to v-face corner (advection of vv)
  RealArr1DHost vx;  // v interpolated to u-face corner (advection of uu)

  // Geometry (areas / volumes).
  RealArr1DHost Vs;    // cell water volume = dept*dx*dy
  RealArr1DHost Vsn;   // previous Vs
  RealArr1DHost Asx;   // x+ face wet area = deptx*dy
  RealArr1DHost Asy;   // y+ face wet area = depty*dx
  RealArr1DHost Asz;   // top area = dx*dy if wet else 0
  RealArr1DHost Aszx;  // 0.5*(Asz + Asz_iP)
  RealArr1DHost Aszy;  // 0.5*(Asz + Asz_jP)
  RealArr1DHost Vsx;   // 0.5*(Vs + Vs_iP)
  RealArr1DHost Vsy;   // 0.5*(Vs + Vs_jP)

  // Drag / momentum.
  RealArr1DHost CDx;   // Manning drag coefficient (x)
  RealArr1DHost CDy;   // Manning drag coefficient (y)
  RealArr1DHost Dx;    // implicit drag factor (x)
  RealArr1DHost Dy;    // implicit drag factor (y)
  RealArr1DHost Ex;    // explicit momentum source (x)
  RealArr1DHost Ey;    // explicit momentum source (y)
  RealArr1DHost cflx;  // x face CFL number
  RealArr1DHost cfly;  // y face CFL number
  RealArr1DHost Fu;    // x face volume flux = uu*Asx
  RealArr1DHost Fv;    // y face volume flux = vv*Asy

  // Forcing.
  RealArr1DHost rain;  // rainfall rate (m/s)
  RealArr1DHost evap;  // evaporation rate (m/s)

  // Linear-system coefficients (legacy Sxp/Sxm/Syp/Sym/Sct/Srhs) and limiter/waterfall.
  RealArr1DHost Sxp;
  RealArr1DHost Sxm;
  RealArr1DHost Syp;
  RealArr1DHost Sym;
  RealArr1DHost Sct;
  RealArr1DHost Srhs;
  RealArr1DHost cfl_active;  // wetting-limiter flag (set by cflLimiter, consumed by velocity)
  RealArr1DHost wtfx;        // waterfall marker (x)
  RealArr1DHost wtfy;        // waterfall marker (y)

  // Every field array, paired with a stable name (used for full-state checkpoint/restart,
  // P7.4). Order is irrelevant; restore matches by name.
  std::vector<std::pair<std::string, RealArr1DHost*>> namedViews() {
    return {{"eta", &eta},     {"etan", &etan},   {"bottom", &bottom},
            {"bottomXP", &bottomXP}, {"bottomYP", &bottomYP}, {"dept", &dept},
            {"deptx", &deptx}, {"depty", &depty}, {"uu", &uu},
            {"vv", &vv},       {"un", &un},       {"vn", &vn},
            {"uy", &uy},       {"vx", &vx},       {"Vs", &Vs},
            {"Vsn", &Vsn},     {"Asx", &Asx},     {"Asy", &Asy},
            {"Asz", &Asz},     {"Aszx", &Aszx},   {"Aszy", &Aszy},
            {"Vsx", &Vsx},     {"Vsy", &Vsy},     {"CDx", &CDx},
            {"CDy", &CDy},     {"Dx", &Dx},       {"Dy", &Dy},
            {"Ex", &Ex},       {"Ey", &Ey},       {"cflx", &cflx},
            {"cfly", &cfly},   {"Fu", &Fu},       {"Fv", &Fv},
            {"rain", &rain},   {"evap", &evap},   {"Sxp", &Sxp},
            {"Sxm", &Sxm},     {"Syp", &Syp},     {"Sym", &Sym},
            {"Sct", &Sct},     {"Srhs", &Srhs},   {"cfl_active", &cfl_active},
            {"wtfx", &wtfx},   {"wtfy", &wtfy}};
  }
};

}  // namespace frehg2

#endif  // FREHG2_SWE_SWE_FIELDS_HPP
