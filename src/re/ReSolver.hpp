#ifndef FREHG2_RE_SOLVER_HPP
#define FREHG2_RE_SOLVER_HPP

#include "core/Grid.hpp"

#include <array>
#include <string>
#include <vector>

namespace frehg2 {

struct SoilParameters {
    real ksx = 0.0;
    real ksy = 0.0;
    real ksz = 0.0;
    real soil_a = 1.0;
    real soil_n = 1.5;
    real wcs = 0.33;
    real wcr = 0.0;
    real aev = -0.02;
};

struct ReParameters {
    real dt = 1.0e-4;
    real dt_min = 1.0e-4;
    real dt_max = 2.0;
    real co_max = 2.0;
    real ksx = 0.0;
    real ksy = 0.0;
    real ksz = 0.0;
    real ss = 0.0;
    real soil_a = 1.0;
    real soil_n = 1.5;
    real wcs = 0.33;
    real wcr = 0.0;
    real init_wc = 0.0;
    real init_h = 1.0;
    real init_wt_rel = 0.0;
    real init_wt_abs = 0.0;
    real qtop = 0.0;
    real qbot = 0.0;
    real htop = 0.0;
    real hbot = 0.0;
    real qyp = 0.0;
    real qym = 0.0;
    real aev = -0.02;
    bool use_vg = true;
    bool use_mvg = false;
    bool use_full3d = false;
    bool use_corrector = true;
    bool dt_adjust = true;
    bool follow_terrain = false;
    std::array<int, 6> bc_type = {0, 0, 0, 0, 0, 1};
    std::vector<real> qtop_surface;
    std::vector<SoilParameters> soil_table;
    std::vector<int> soil_id;
};

struct ReMatrixEntry {
    int row = 0;
    int col = 0;
    real value = 0.0;
};

struct ReLinearSystem {
    int n = 0;
    std::vector<ReMatrixEntry> entries;
    std::vector<real> rhs;
};

struct LegacyReState {
    std::vector<real> h;
    std::vector<real> hn;
    std::vector<real> hnm;
    std::vector<real> wc;
    std::vector<real> wcn;
    std::vector<real> wch;
    std::vector<real> hwc;
    std::vector<real> ch;
    std::vector<real> kx;
    std::vector<real> ky;
    std::vector<real> kz;
    std::vector<real> qx;
    std::vector<real> qy;
    std::vector<real> qz;
    real qz_top = 0.0;
    std::vector<real> qz_top_surface;
    std::vector<real> gxp;
    std::vector<real> gxm;
    std::vector<real> gyp;
    std::vector<real> gym;
    std::vector<real> gzp;
    std::vector<real> gzm;
    std::vector<real> gct;
    std::vector<real> grhs;
    std::vector<real> bot3d;
    std::vector<real> dz3d;
    std::vector<real> zcntr;
    std::vector<int> active;
    std::vector<int> istop;
    std::vector<int> soil_id;
    real dtg = 1.0e-4;
};

class ReSolver {
public:
    ReSolver(Grid grid, ReParameters parameters);

    const Grid& grid() const noexcept;
    const ReParameters& parameters() const noexcept;

    static int legacyNz(real bath, real bot_z, real dz, real dz_multiplier);
    static real waterContentFromHead(real h, const ReParameters& parameters);
    static real headFromWaterContent(real wc, const ReParameters& parameters);
    static real specificCapacity(real h, const ReParameters& parameters);
    static real hydraulicConductivity(real h, real ks, const ReParameters& parameters);
    static real conductivityFace(real h_left, real h_right, real ks, const ReParameters& parameters);
    static real conductivityDerivativeWc(real wc, real ks, const ReParameters& parameters);
    static SoilParameters uniformSoil(const ReParameters& parameters);
    static ReParameters parametersForSoil(const ReParameters& parameters, const SoilParameters& soil);
    static std::vector<int> expandSurfaceSoilMap(
        const Grid& grid,
        const std::vector<int>& surface_soil_id);
    static std::vector<int> readSoilMap(const Grid& grid, const std::string& filename);
    static std::vector<int> readSurfaceSoilMap(const Grid& grid, const std::string& filename);

    LegacyReState initializeLegacyState(real bath, real bot_z) const;
    void advanceLegacyStep(LegacyReState& state) const;
    void computeConductivityFaces(LegacyReState& state) const;
    ReLinearSystem assemblePredictorSystem(LegacyReState& state) const;
    std::vector<real> solveLinearSystem(const ReLinearSystem& system) const;
    void computeDarcyFlux(LegacyReState& state) const;
    void updateWaterContent(LegacyReState& state) const;
    void enforceLegacyMoistureConsistency(LegacyReState& state) const;
    real adaptTimeStep(const LegacyReState& state) const;

private:
    Grid grid_;
    ReParameters parameters_;

    int index(int i, int j, int k) const;
    int zp(int cell) const;
    int zm(int cell) const;
    real qtopForSurfaceCell(int surface_cell) const;
    const SoilParameters& soilForCell(const LegacyReState& state, int cell) const;
    ReParameters parametersForCell(const LegacyReState& state, int cell) const;
    bool isTop(int cell) const;
    bool isAdjacentToSaturation(const LegacyReState& state, int cell) const;
};

}  // namespace frehg2

#endif  // FREHG2_RE_SOLVER_HPP
