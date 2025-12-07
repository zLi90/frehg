#ifndef FREHG_STATE_VARIABLES_HPP
#define FREHG_STATE_VARIABLES_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "ActiveCellMesh.hpp"
#include <memory>

// ============================================================================
//                      BASE STATE VARIABLES CLASS
// ============================================================================
// Shared state variables for both surface water and groundwater domains

class BaseStateVariables {
public:
    // --- Pressure/Head Variables ---
    View1D<Scalar> pressure;        // Pressure/head (eta for SW, h for GW)
    View1D<Scalar> pressure_old;     // Previous time step pressure
    
    // --- Velocity Variables (shared) ---
    View1D<Scalar> velocity_x;      // Velocity in x-direction (u)
    View1D<Scalar> velocity_y;       // Velocity in y-direction (v)
    View1D<Scalar> velocity_x_old;  // Previous time step velocity x
    View1D<Scalar> velocity_y_old;  // Previous time step velocity y
    
    // --- Flux Variables ---
    View1D<Scalar> flux_x;           // Flux in x-direction
    View1D<Scalar> flux_y;           // Flux in y-direction
    
    // --- Matrix Coefficients (for linear systems) ---
    View1D<Scalar> matrix_diag;      // Diagonal coefficient
    View1D<Scalar> matrix_rhs;       // Right-hand side
    
    // Constructor
    BaseStateVariables(Ordinal num_cells) {
        pressure = View1D<Scalar>("pressure", num_cells);
        pressure_old = View1D<Scalar>("pressure_old", num_cells);
        velocity_x = View1D<Scalar>("velocity_x", num_cells);
        velocity_y = View1D<Scalar>("velocity_y", num_cells);
        velocity_x_old = View1D<Scalar>("velocity_x_old", num_cells);
        velocity_y_old = View1D<Scalar>("velocity_y_old", num_cells);
        flux_x = View1D<Scalar>("flux_x", num_cells);
        flux_y = View1D<Scalar>("flux_y", num_cells);
        matrix_diag = View1D<Scalar>("matrix_diag", num_cells);
        matrix_rhs = View1D<Scalar>("matrix_rhs", num_cells);
    }
    
    virtual ~BaseStateVariables() = default;
    
    // Update old values (called at start of time step)
    virtual void update_old_values() {
        Kokkos::deep_copy(pressure_old, pressure);
        Kokkos::deep_copy(velocity_x_old, velocity_x);
        Kokkos::deep_copy(velocity_y_old, velocity_y);
    }
};

// ============================================================================
//                  SURFACE WATER STATE VARIABLES
// ============================================================================

class SwStateVariables : public BaseStateVariables {
public:
    // --- Surface Elevation (pressure = eta) ---
    // pressure is eta (surface elevation)
    // pressure_old is etan (previous surface elevation)
    
    // --- Depth Variables ---
    View1D<Scalar> depth;            // Total water depth (eta - bottom)
    View1D<Scalar> depth_x;          // Depth at x-face (deptx)
    View1D<Scalar> depth_y;          // Depth at y-face (depty)
    
    // --- Bathymetry ---
    View1D<Scalar> bottom;           // Bathymetry/elevation (z)
    View1D<Scalar> bottom_x;         // Bathymetry at x+ face (bottomXP)
    View1D<Scalar> bottom_y;         // Bathymetry at y+ face (bottomYP)
    
    // --- Interpolated Velocities ---
    View1D<Scalar> velocity_x_at_y;  // u interpolated to y-face (uy)
    View1D<Scalar> velocity_y_at_x;  // v interpolated to x-face (vx)
    
    // --- Momentum Source Terms ---
    View1D<Scalar> momentum_x;       // Explicit momentum term X (Ex)
    View1D<Scalar> momentum_y;       // Explicit momentum term Y (Ey)
    
    // --- Drag Coefficients ---
    View1D<Scalar> drag_factor_x;    // Drag factor X (Dx)
    View1D<Scalar> drag_factor_y;    // Drag factor Y (Dy)
    View1D<Scalar> drag_coef_x;      // Drag coefficient X (CDx)
    View1D<Scalar> drag_coef_y;      // Drag coefficient Y (CDy)
    
    // --- Geometry (Face Areas & Volumes) ---
    View1D<Scalar> area_x;           // Cross-sectional area X-face (Asx)
    View1D<Scalar> area_y;           // Cross-sectional area Y-face (Asy)
    View1D<Scalar> area_top;         // Top area / cell area (Asz)
    View1D<Scalar> area_top_x;       // Top area at x-face (Aszx)
    View1D<Scalar> area_top_y;       // Top area at y-face (Aszy)
    
    View1D<Scalar> volume;           // Cell volume (Vs)
    View1D<Scalar> volume_old;       // Previous volume (Vsn)
    View1D<Scalar> volume_flux;      // Volume from flux (Vflux)
    View1D<Scalar> volume_x;         // Staggered volume X (Vsx)
    View1D<Scalar> volume_y;         // Staggered volume Y (Vsy)
    
    // --- Matrix Coefficients (Surface Water) ---
    View1D<Scalar> matrix_xp;        // X+ coefficient (Sxp)
    View1D<Scalar> matrix_xm;        // X- coefficient (Sxm)
    View1D<Scalar> matrix_yp;        // Y+ coefficient (Syp)
    View1D<Scalar> matrix_ym;        // Y- coefficient (Sym)
    
    // --- Subsurface Coupling ---
    View1D<Scalar> seepage;          // Seepage flux (qseepage)
    View1D<Scalar> seepage_old;      // Previous seepage (qseepage_old)
    View1D<Scalar> seepage_rate;     // Seepage rate (qss)
    View1D<int> reset_seepage;       // Flag to reset seepage
    
    // --- Stability & Limiting ---
    View1D<Scalar> cfl_x;            // CFL number in x (cflx)
    View1D<Scalar> cfl_y;            // CFL number in y (cfly)
    View1D<int> cfl_active;          // CFL limiter active flag
    
    // --- Waterfall Detection ---
    View1D<int> waterfall_x;        // Waterfall flag x (wtfx)
    View1D<int> waterfall_y;        // Waterfall flag y (wtfy)
    
    // Constructor
    SwStateVariables(Ordinal num_cells) : BaseStateVariables(num_cells) {
        // Allocate all surface water specific variables
        depth = View1D<Scalar>("depth", num_cells);
        depth_x = View1D<Scalar>("depth_x", num_cells);
        depth_y = View1D<Scalar>("depth_y", num_cells);
        
        bottom = View1D<Scalar>("bottom", num_cells);
        bottom_x = View1D<Scalar>("bottom_x", num_cells);
        bottom_y = View1D<Scalar>("bottom_y", num_cells);
        
        velocity_x_at_y = View1D<Scalar>("velocity_x_at_y", num_cells);
        velocity_y_at_x = View1D<Scalar>("velocity_y_at_x", num_cells);
        
        momentum_x = View1D<Scalar>("momentum_x", num_cells);
        momentum_y = View1D<Scalar>("momentum_y", num_cells);
        
        drag_factor_x = View1D<Scalar>("drag_factor_x", num_cells);
        drag_factor_y = View1D<Scalar>("drag_factor_y", num_cells);
        drag_coef_x = View1D<Scalar>("drag_coef_x", num_cells);
        drag_coef_y = View1D<Scalar>("drag_coef_y", num_cells);
        
        area_x = View1D<Scalar>("area_x", num_cells);
        area_y = View1D<Scalar>("area_y", num_cells);
        area_top = View1D<Scalar>("area_top", num_cells);
        area_top_x = View1D<Scalar>("area_top_x", num_cells);
        area_top_y = View1D<Scalar>("area_top_y", num_cells);
        
        volume = View1D<Scalar>("volume", num_cells);
        volume_old = View1D<Scalar>("volume_old", num_cells);
        volume_flux = View1D<Scalar>("volume_flux", num_cells);
        volume_x = View1D<Scalar>("volume_x", num_cells);
        volume_y = View1D<Scalar>("volume_y", num_cells);
        
        matrix_xp = View1D<Scalar>("matrix_xp", num_cells);
        matrix_xm = View1D<Scalar>("matrix_xm", num_cells);
        matrix_yp = View1D<Scalar>("matrix_yp", num_cells);
        matrix_ym = View1D<Scalar>("matrix_ym", num_cells);
        
        seepage = View1D<Scalar>("seepage", num_cells);
        seepage_old = View1D<Scalar>("seepage_old", num_cells);
        seepage_rate = View1D<Scalar>("seepage_rate", num_cells);
        reset_seepage = View1D<int>("reset_seepage", num_cells);
        
        cfl_x = View1D<Scalar>("cfl_x", num_cells);
        cfl_y = View1D<Scalar>("cfl_y", num_cells);
        cfl_active = View1D<int>("cfl_active", num_cells);
        
        waterfall_x = View1D<int>("waterfall_x", num_cells);
        waterfall_y = View1D<int>("waterfall_y", num_cells);
    }
    
    // Override update_old_values to include surface water specific
    void update_old_values() override {
        BaseStateVariables::update_old_values();
        Kokkos::deep_copy(volume_old, volume);
        Kokkos::deep_copy(seepage_old, seepage);
    }
};

// ============================================================================
//                  GROUNDWATER STATE VARIABLES
// ============================================================================

class GwStateVariables : public BaseStateVariables {
public:
    // --- Head Variables (pressure = h) ---
    // pressure is h (hydraulic head)
    // pressure_old is hn (previous head)
    View1D<Scalar> head_prev_prev;  // Previous-previous head (hnm)
    View1D<Scalar> head_predictor;  // Predictor head (hp)
    View1D<Scalar> head_increment;  // Head increment for Newton (h_incr)
    
    // --- Water Content Variables ---
    View1D<Scalar> water_content;    // Water content (wc)
    View1D<Scalar> water_content_old;// Previous water content (wcn)
    View1D<Scalar> water_content_from_head; // wc from h (wch)
    View1D<Scalar> head_from_water_content; // h from wc (hwc)
    
    // --- Saturated Properties ---
    View1D<Scalar> water_content_sat; // Saturated water content (wcs)
    View1D<Scalar> water_content_res; // Residual water content (wcr)
    
    // --- Van Genuchten Parameters ---
    View1D<Scalar> vg_alpha;          // Van Genuchten alpha parameter [1/L]
    View1D<Scalar> vg_n;              // Van Genuchten n parameter (dimensionless)
    View1D<Scalar> vg_m;              // Van Genuchten m parameter (m = 1 - 1/n)
    View1D<Scalar> vg_ha;             // Air entry value [L] (ha >= 0, if ha=0 uses original vG model)
    
    // --- Hydraulic Conductivity ---
    View1D<Scalar> conductivity_x;  // Kx (face conductivity)
    View1D<Scalar> conductivity_y;  // Ky (face conductivity)
    View1D<Scalar> conductivity_z;  // Kz (face conductivity)
    View1D<Scalar> conductivity_sat_x; // Saturated Kx (Ksx)
    View1D<Scalar> conductivity_sat_y; // Saturated Ky (Ksy)
    View1D<Scalar> conductivity_sat_z; // Saturated Kz (Ksz)
    
    // --- Flux Variables (3D) ---
    View1D<Scalar> flux_z;           // Flux in z-direction (qz)
    
    // --- Specific Storage ---
    View1D<Scalar> specific_storage; // Specific storage coefficient (ch)
    
    // --- Matrix Coefficients (Groundwater) ---
    View1D<Scalar> matrix_xp;        // X+ coefficient (Gxp)
    View1D<Scalar> matrix_xm;        // X- coefficient (Gxm)
    View1D<Scalar> matrix_yp;        // Y+ coefficient (Gyp)
    View1D<Scalar> matrix_ym;        // Y- coefficient (Gym)
    View1D<Scalar> matrix_zp;        // Z+ coefficient (Gzp)
    View1D<Scalar> matrix_zm;        // Z- coefficient (Gzm)
    
    // --- Baroclinic Variables (density/viscosity) ---
    View1D<Scalar> density_ratio;    // Density ratio (r_rho)
    View1D<Scalar> density_ratio_old; // Previous density ratio (r_rhon)
    View1D<Scalar> viscosity_ratio; // Viscosity ratio (r_visc)
    
    // Face density/viscosity ratios
    View1D<Scalar> density_ratio_xp; // r_rhoxp
    View1D<Scalar> density_ratio_yp; // r_rhoyp
    View1D<Scalar> density_ratio_zp; // r_rhozp
    View1D<Scalar> viscosity_ratio_xp; // r_viscxp
    View1D<Scalar> viscosity_ratio_yp; // r_viscyp
    View1D<Scalar> viscosity_ratio_zp; // r_visczp
    
    // --- Volume Variables ---
    View1D<Scalar> volume;           // Cell volume (Vg)
    View1D<Scalar> volume_old;       // Previous volume (Vgn)
    View1D<Scalar> volume_flux;      // Volume from flux (Vgflux)
    View1D<Scalar> volume_loss;      // Volume loss (vloss)
    View1D<Scalar> room;             // Available room for water (room)
    
    // --- Seepage/Coupling ---
    View1D<Scalar> seepage_top;      // Top boundary flux (qtop)
    View1D<Scalar> seepage_bot;      // Bottom boundary flux (qbot)
    
    // --- Newton Iteration Variables ---
    View1D<Scalar> residual;         // Residual (resi)
    View1D<Scalar> head_gradient;    // Head gradient array (dh6, size 6)
    View1D<Scalar> moisture_split;   // Moisture split ratio (rsplit, size 6)
    View1D<int> repeat;              // Repeat flag for allocation
    
    // Constructor
    GwStateVariables(Ordinal num_cells) : BaseStateVariables(num_cells) {
        // Allocate all groundwater specific variables
        head_prev_prev = View1D<Scalar>("head_prev_prev", num_cells);
        head_predictor = View1D<Scalar>("head_predictor", num_cells);
        head_increment = View1D<Scalar>("head_increment", num_cells);
        
        water_content = View1D<Scalar>("water_content", num_cells);
        water_content_old = View1D<Scalar>("water_content_old", num_cells);
        water_content_from_head = View1D<Scalar>("water_content_from_head", num_cells);
        head_from_water_content = View1D<Scalar>("head_from_water_content", num_cells);
        
        water_content_sat = View1D<Scalar>("water_content_sat", num_cells);
        water_content_res = View1D<Scalar>("water_content_res", num_cells);
        
        // Van Genuchten parameters (default values: typical for sandy loam)
        vg_alpha = View1D<Scalar>("vg_alpha", num_cells);
        vg_n = View1D<Scalar>("vg_n", num_cells);
        vg_m = View1D<Scalar>("vg_m", num_cells);
        vg_ha = View1D<Scalar>("vg_ha", num_cells);
        
        conductivity_x = View1D<Scalar>("conductivity_x", num_cells);
        conductivity_y = View1D<Scalar>("conductivity_y", num_cells);
        conductivity_z = View1D<Scalar>("conductivity_z", num_cells);
        conductivity_sat_x = View1D<Scalar>("conductivity_sat_x", num_cells);
        conductivity_sat_y = View1D<Scalar>("conductivity_sat_y", num_cells);
        conductivity_sat_z = View1D<Scalar>("conductivity_sat_z", num_cells);
        
        flux_z = View1D<Scalar>("flux_z", num_cells);
        
        specific_storage = View1D<Scalar>("specific_storage", num_cells);
        
        matrix_xp = View1D<Scalar>("matrix_xp", num_cells);
        matrix_xm = View1D<Scalar>("matrix_xm", num_cells);
        matrix_yp = View1D<Scalar>("matrix_yp", num_cells);
        matrix_ym = View1D<Scalar>("matrix_ym", num_cells);
        matrix_zp = View1D<Scalar>("matrix_zp", num_cells);
        matrix_zm = View1D<Scalar>("matrix_zm", num_cells);
        
        density_ratio = View1D<Scalar>("density_ratio", num_cells);
        density_ratio_old = View1D<Scalar>("density_ratio_old", num_cells);
        viscosity_ratio = View1D<Scalar>("viscosity_ratio", num_cells);
        
        density_ratio_xp = View1D<Scalar>("density_ratio_xp", num_cells);
        density_ratio_yp = View1D<Scalar>("density_ratio_yp", num_cells);
        density_ratio_zp = View1D<Scalar>("density_ratio_zp", num_cells);
        viscosity_ratio_xp = View1D<Scalar>("viscosity_ratio_xp", num_cells);
        viscosity_ratio_yp = View1D<Scalar>("viscosity_ratio_yp", num_cells);
        viscosity_ratio_zp = View1D<Scalar>("viscosity_ratio_zp", num_cells);
        
        volume = View1D<Scalar>("volume", num_cells);
        volume_old = View1D<Scalar>("volume_old", num_cells);
        volume_flux = View1D<Scalar>("volume_flux", num_cells);
        volume_loss = View1D<Scalar>("volume_loss", num_cells);
        room = View1D<Scalar>("room", num_cells);
        
        seepage_top = View1D<Scalar>("seepage_top", num_cells);
        seepage_bot = View1D<Scalar>("seepage_bot", num_cells);
        
        residual = View1D<Scalar>("residual", num_cells);
        head_gradient = View1D<Scalar>("head_gradient", 6 * num_cells); // 6 directions
        moisture_split = View1D<Scalar>("moisture_split", 6 * num_cells); // 6 directions
        repeat = View1D<int>("repeat", num_cells);
        
        // Initialize van Genuchten parameters with default values
        // (can be overridden during initialization from input files)
        // Default values: typical for sandy loam
        auto _vg_alpha = vg_alpha;
        auto _vg_n = vg_n;
        auto _vg_m = vg_m;
        auto _vg_ha = vg_ha;
        using RangePolicy = Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>;
        Kokkos::parallel_for(RangePolicy(0, num_cells),
            KOKKOS_LAMBDA (const Ordinal i) {
                _vg_alpha(i) = 0.01;  // Default: 0.01 cm^-1 (typical for sandy loam)
                _vg_n(i) = 1.5;       // Default: n = 1.5
                _vg_m(i) = 1.0 - 1.0 / 1.5;  // m = 1 - 1/n ≈ 0.333
                _vg_ha(i) = 0.0;     // Default: 0.0 (original van Genuchten model)
            });
    }
    
    // Override update_old_values to include groundwater specific
    void update_old_values() override {
        BaseStateVariables::update_old_values();
        Kokkos::deep_copy(head_prev_prev, pressure_old);
        Kokkos::deep_copy(water_content_old, water_content);
        Kokkos::deep_copy(volume_old, volume);
        Kokkos::deep_copy(density_ratio_old, density_ratio);
    }
    
    // Access head gradient for direction dir (0-5: xm, xp, ym, yp, zm, zp)
    KOKKOS_INLINE_FUNCTION
    Scalar& head_grad(Ordinal cell_idx, Ordinal dir) {
        return head_gradient(cell_idx * 6 + dir);
    }
    
    // Access moisture split for direction dir
    KOKKOS_INLINE_FUNCTION
    Scalar& moisture_split_dir(Ordinal cell_idx, Ordinal dir) {
        return moisture_split(cell_idx * 6 + dir);
    }
};

#endif // FREHG_STATE_VARIABLES_HPP

