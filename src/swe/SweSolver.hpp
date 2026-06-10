#ifndef FREHG2_SWE_SOLVER_HPP
#define FREHG2_SWE_SOLVER_HPP

#include "core/Domain.hpp"
#include "core/State.hpp"

#include <vector>

namespace frehg2 {

struct SweFreeOutflowBoundary {
    std::vector<index_t> cells;
    real normal_x = 1.0;
    real normal_y = 0.0;
};

struct SweParameters {
    real dt = 1.0;
    real gravity = 9.81;
    real min_depth = 1.0e-8;
    real manning = 0.0;
    real rain_rate = 0.0;
    real evap_rate = 0.0;
    real cfl = 1.0;
    real viscx = 1.0e-6;
    real viscy = 1.0e-6;
    real wtfh = 1.0e-8;
    real hD = 0.1;
    std::vector<SweFreeOutflowBoundary> free_outflow_boundaries;
};

struct SweMatrixEntry {
    int row = 0;
    int col = 0;
    real value = 0.0;
};

struct SweLinearSystem {
    int n = 0;
    std::vector<SweMatrixEntry> entries;
    std::vector<real> rhs;
    std::vector<int> dom2mat;
    std::vector<int> mat2dom;
};

struct LegacySweState {
    std::vector<real> bottom;
    std::vector<real> bottom_xp;
    std::vector<real> bottom_yp;
    std::vector<int> active;

    std::vector<real> eta;
    std::vector<real> etan;
    std::vector<real> dept;
    std::vector<real> deptx;
    std::vector<real> depty;

    std::vector<real> uu;
    std::vector<real> un;
    std::vector<real> uy;
    std::vector<real> vv;
    std::vector<real> vn;
    std::vector<real> vx;

    std::vector<real> Fu;
    std::vector<real> Fv;
    std::vector<real> Ex;
    std::vector<real> Ey;
    std::vector<real> Dx;
    std::vector<real> Dy;
    std::vector<real> CDx;
    std::vector<real> CDy;

    std::vector<real> Vs;
    std::vector<real> Vsn;
    std::vector<real> Vflux;
    std::vector<real> Vsx;
    std::vector<real> Vsy;
    std::vector<real> Asz;
    std::vector<real> Aszx;
    std::vector<real> Aszy;
    std::vector<real> Asx;
    std::vector<real> Asy;

    std::vector<real> Sxp;
    std::vector<real> Sxm;
    std::vector<real> Syp;
    std::vector<real> Sym;
    std::vector<real> Sct;
    std::vector<real> Srhs;

    std::vector<real> wtfx;
    std::vector<real> wtfy;
    std::vector<real> cflx;
    std::vector<real> cfly;
    std::vector<real> cfl_active;
    real rain_sum = 0.0;
};

class SweSolver {
public:
    SweSolver(Grid grid, SweParameters parameters);

    const Grid& grid() const noexcept;
    const SweParameters& parameters() const noexcept;

    void initialize(const Domain& domain, State& state, real initial_eta) const;

    SweLinearSystem assembleLinearSystem(
        const std::vector<real>& eta,
        const std::vector<real>& bed,
        const std::vector<int>& active_mask,
        const std::vector<real>& ex,
        const std::vector<real>& ey,
        const std::vector<real>& dx_factor,
        const std::vector<real>& dy_factor) const;

    std::vector<real> solveLinearSystem(const SweLinearSystem& system) const;

    void applyRainEvap(std::vector<real>& eta, const std::vector<real>& bed) const;
    real computeCflDiagnostic(
        const std::vector<real>& h,
        const std::vector<real>& u,
        const std::vector<real>& v) const;
    void updateVelocity(
        std::vector<real>& u,
        std::vector<real>& v,
        const std::vector<real>& eta_old,
        const std::vector<real>& eta_new) const;

    LegacySweState initializeLegacyState(const std::vector<real>& bed, real initial_eta) const;
    void advanceLegacyStep(LegacySweState& state, real rain_rate, real evap_rate) const;

private:
    Grid grid_;
    SweParameters parameters_;

    int xp(int i, int j) const;
    int xm(int i, int j) const;
    int yp(int i, int j) const;
    int ym(int i, int j) const;

    void enforceLegacySurfaceBc(LegacySweState& state) const;
    void enforceLegacyVelocityBc(LegacySweState& state) const;
    void updateLegacyDepth(LegacySweState& state) const;
    void updateLegacySubgridVariables(LegacySweState& state, bool zero_degenerate_faces) const;
    void updateLegacyDragCoefficient(LegacySweState& state) const;
    void legacyMomentumSource(LegacySweState& state) const;
    bool isFreeOutflowBoundaryCell(index_t cell, real normal_x, real normal_y) const;
    void enforceLegacyReflectiveBoundaryFluxes(LegacySweState& state) const;
    SweLinearSystem assembleLegacyLinearSystem(LegacySweState& state) const;
    void legacyCflLimiter(LegacySweState& state, real rain_rate) const;
    void legacyEvapRain(LegacySweState& state, real rain_rate, real evap_rate) const;
    void legacyWaterfallLocation(LegacySweState& state) const;
    void updateLegacyVelocity(LegacySweState& state) const;
    void applyLegacyFreeOutflowBc(LegacySweState& state) const;
    void interpolateLegacyVelocity(LegacySweState& state) const;
    void fillLegacyGhosts(LegacySweState& state) const;
};

}  // namespace frehg2

#endif  // FREHG2_SWE_SOLVER_HPP
