#ifndef FREHG_GROUNDWATER_SOLVER_HPP
#define FREHG_GROUNDWATER_SOLVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "PCGSolver.hpp"
#include "BoundaryConditions.hpp"
#include <memory>
#include <cmath>

// ============================================================================
//                      GROUNDWATER EQUATION SOLVER
// ============================================================================
// Solves the 3D Richards equation for variably saturated groundwater flow
// Uses Kokkos for parallelization and PCG for linear system solution
// Implements predictor-corrector scheme (no Newton iteration)

class GroundwaterSolver {
public:
    // --- Domain and Mesh ---
    GwDomain& domain;
    ActiveCellMesh& active_mesh;
    
    // --- State Variables ---
    GwStateVariables& state;
    
    // --- Physical Parameters ---
    Scalar dt;                  // Time step
    Scalar Ss;                  // Specific storage
    Scalar min_depth;           // Minimum depth threshold
    bool use_corrector;         // Use corrector step
    bool follow_terrain;        // Follow terrain (slope effects)
    bool baroclinic;            // Baroclinic effects (density/viscosity)
    
    // --- Linear Solver ---
    std::unique_ptr<PCGSolver3D> pcg_solver;
    PreconditionerType precond_type;
    Scalar solver_tolerance;
    Ordinal max_solver_iterations;
    
    // --- Solution Vector (active cells only) ---
    View1D<Scalar> solution;    // Solution vector for PCG (h_new in active space)
    
    // --- Boundary Conditions ---
    GwBoundaryConditionManager* bc_manager_;  // Pointer to boundary condition manager
    
    // Constructor
    GroundwaterSolver(GwDomain& _domain,
                      ActiveCellMesh& _active_mesh,
                      GwStateVariables& _state,
                      Scalar _dt,
                      Scalar _Ss = 1.0e-5,
                      Scalar _min_depth = 1.0e-6,
                      bool _use_corrector = true,
                      bool _follow_terrain = false,
                      bool _baroclinic = false,
                      PreconditionerType _precond_type = PreconditionerType::JACOBI,
                      Scalar _solver_tolerance = 1e-8,
                      Ordinal _max_solver_iterations = 10000,
                      GwBoundaryConditionManager* _bc_manager = nullptr)
        : domain(_domain), active_mesh(_active_mesh), state(_state),
          dt(_dt), Ss(_Ss), min_depth(_min_depth),
          use_corrector(_use_corrector), follow_terrain(_follow_terrain),
          baroclinic(_baroclinic),
          precond_type(_precond_type), solver_tolerance(_solver_tolerance),
          max_solver_iterations(_max_solver_iterations),
          bc_manager_(_bc_manager) {
        
        // Initialize PCG solver
        pcg_solver = std::make_unique<PCGSolver3D>(
            active_mesh.num_active, precond_type, solver_tolerance, max_solver_iterations);
        
        // Allocate solution vector for active cells
        solution = View1D<Scalar>("solution", active_mesh.num_active);
    }
    
    // Set boundary condition manager
    void set_boundary_condition_manager(GwBoundaryConditionManager* bc_manager) {
        bc_manager_ = bc_manager;
    }

    // Main solver function (predictor-corrector scheme)
    void solve(Scalar current_time = 0.0) {
        // Save old values
        state.update_old_values(domain.num_cells_3d_total);
        
        // Compute water content from head and specific moisture capacity
        compute_water_content_from_head();
        compute_specific_moisture_capacity();
        
        // >>> PREDICTOR STEP <<<
        // 1. Compute hydraulic conductivity on faces
        compute_K_face();
        
        // 2. Compute baroclinic effects (density/viscosity ratios)
        if (baroclinic) {
            baroclinic_face();
        }
        
        // 3. Assemble linear system
        assemble_linear_system();
        
        // 4. Apply boundary conditions to matrix
        apply_boundary_conditions_to_matrix(current_time);
        
        // 5. Solve for new head using PCG
        solve_linear_system();
        
        // 6. Enforce boundary conditions on solution
        enforce_boundary_conditions_on_solution(current_time);
        
        // >>> CORRECTOR STEP <<<
        if (use_corrector) {
            // 1. Recompute hydraulic conductivity
            compute_K_face();
            
            // 2. Compute groundwater fluxes
            compute_groundwater_flux();
            
            // 3. Check available room (pore space)
            check_room();
            
            // 4. Update water content from fluxes
            update_water_content();
        } else {
            // Direct update from head
            compute_water_content_from_head();
        }
        
        // Final update of derived quantities
        compute_head_from_water_content();
        compute_water_content_from_head();
    }

private:
    // ========================================================================
    // HELPER FUNCTIONS (Soil Property Functions - Placeholders)
    // ========================================================================
    // These functions should be implemented based on soil property models
    // For now, we assume they are available or will be implemented later
    
    // Compute hydraulic conductivity from head and saturated conductivity
    // K = Ks * Kr(h), where Kr is relative permeability
    KOKKOS_INLINE_FUNCTION
    Scalar compute_K(Scalar h, Scalar Ks, Ordinal idx) const {
        // Placeholder: assume K = Ks for saturated, reduce for unsaturated
        // This should be replaced with actual soil property model
        if (h >= 0.0) {
            return Ks;
        } else {
            // Simplified: exponential reduction for unsaturated
            Scalar wc = state.water_content(idx);
            Scalar wcs = state.water_content_sat(idx);
            Scalar wcr = state.water_content_res(idx);
            if (wcs > wcr) {
                Scalar Se = (wc - wcr) / (wcs - wcr); // Effective saturation
                return Ks * Se * Se * Se; // Cubic relationship (simplified)
            }
            return Ks * 1.0e-6; // Very small for dry
        }
    }
    
    // Compute water content from head using water retention curve
    KOKKOS_INLINE_FUNCTION
    Scalar compute_wch(Scalar h, Ordinal idx) const {
        // Placeholder: simplified water retention curve
        // This should be replaced with van Genuchten or other model
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        
        if (h >= 0.0) {
            return wcs;
        } else {
            // Simplified exponential relationship
            Scalar Se = std::exp(0.1634 * h);
            return wcr + (wcs - wcr) * Se;
        }
    }
    
    // Compute head from water content
    KOKKOS_INLINE_FUNCTION
    Scalar compute_hwc(Ordinal idx) const {
        // Placeholder: inverse of water retention curve
        Scalar wc = state.water_content(idx);
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        
        if (wc >= wcs) {
            return 0.0; // Saturated
        } else if (wc <= wcr) {
            return -1000.0; // Very dry
        } else {
            Scalar Se = (wc - wcr) / (wcs - wcr);
            return std::log(Se) / 0.1634;
        }
    }
    
    // Compute specific moisture capacity (dwc/dh)
    KOKKOS_INLINE_FUNCTION
    Scalar compute_ch(Ordinal idx) const {
        // Placeholder: derivative of water retention curve
        Scalar h = state.pressure(idx);
        Scalar wcs = state.water_content_sat(idx);
        Scalar wcr = state.water_content_res(idx);
        
        if (h >= 0.0) {
            return 0.0; // Saturated
        } else {
            Scalar Se = std::exp(0.1634 * h);
            return 0.1634 * (wcs - wcr) * Se;
        }
    }
    
    // ========================================================================
    // COMPUTE WATER CONTENT FROM HEAD
    // ========================================================================
    void compute_water_content_from_head() {
        auto _pressure = state.pressure;
        auto _water_content = state.water_content;
        auto _water_content_from_head = state.water_content_from_head;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) > 0) {
                    Scalar h = _pressure(domain_idx);
                    _water_content_from_head(domain_idx) = compute_wch(h, domain_idx);
                }
            });
    }
    
    // ========================================================================
    // COMPUTE HEAD FROM WATER CONTENT
    // ========================================================================
    void compute_head_from_water_content() {
        auto _head_from_water_content = state.head_from_water_content;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) > 0) {
                    _head_from_water_content(domain_idx) = compute_hwc(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // COMPUTE SPECIFIC MOISTURE CAPACITY
    // ========================================================================
    void compute_specific_moisture_capacity() {
        // This is computed on-the-fly in matrix assembly
        // Placeholder function for now
    }
    
    // ========================================================================
    // COMPUTE HYDRAULIC CONDUCTIVITY ON FACES
    // ========================================================================
    void compute_K_face() {
        auto _pressure = state.pressure;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _conductivity_sat_x = state.conductivity_sat_x;
        auto _conductivity_sat_y = state.conductivity_sat_y;
        auto _conductivity_sat_z = state.conductivity_sat_z;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _conductivity_x(domain_idx) = 0.0;
                    _conductivity_y(domain_idx) = 0.0;
                    _conductivity_z(domain_idx) = 0.0;
                    return;
                }
                
                // X-direction (face between i and i+1)
                Ordinal active_right = _neighbor_right(i);
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    Scalar h_here = _pressure(domain_idx);
                    Scalar h_right = _pressure(idx_right);
                    Scalar Ks_here = _conductivity_sat_x(domain_idx);
                    Scalar Ks_right = _conductivity_sat_x(idx_right);
                    Scalar K_here = compute_K(h_here, Ks_here, domain_idx);
                    Scalar K_right = compute_K(h_right, Ks_right, idx_right);
                    _conductivity_x(domain_idx) = 0.5 * (K_here + K_right);
                } else {
                    _conductivity_x(domain_idx) = compute_K(_pressure(domain_idx), 
                                                             _conductivity_sat_x(domain_idx), domain_idx);
                }
                
                // Y-direction (face between j and j+1)
                Ordinal active_front = _neighbor_front(i);
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    Scalar h_here = _pressure(domain_idx);
                    Scalar h_front = _pressure(idx_front);
                    Scalar Ks_here = _conductivity_sat_y(domain_idx);
                    Scalar Ks_front = _conductivity_sat_y(idx_front);
                    Scalar K_here = compute_K(h_here, Ks_here, domain_idx);
                    Scalar K_front = compute_K(h_front, Ks_front, idx_front);
                    _conductivity_y(domain_idx) = 0.5 * (K_here + K_front);
                } else {
                    _conductivity_y(domain_idx) = compute_K(_pressure(domain_idx), 
                                                             _conductivity_sat_y(domain_idx), domain_idx);
                }
                
                // Z-direction (face between k and k+1)
                Ordinal active_top = _neighbor_top(i);
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    Scalar h_here = _pressure(domain_idx);
                    Scalar h_top = _pressure(idx_top);
                    Scalar Ks_here = _conductivity_sat_z(domain_idx);
                    Scalar Ks_top = _conductivity_sat_z(idx_top);
                    Scalar K_here = compute_K(h_here, Ks_here, domain_idx);
                    Scalar K_top = compute_K(h_top, Ks_top, idx_top);
                    _conductivity_z(domain_idx) = 0.5 * (K_here + K_top);
                } else {
                    _conductivity_z(domain_idx) = compute_K(_pressure(domain_idx), 
                                                             _conductivity_sat_z(domain_idx), domain_idx);
                }
            });
    }
    
    // ========================================================================
    // COMPUTE BAROCLINIC EFFECTS (Density/Viscosity Ratios)
    // ========================================================================
    void baroclinic_face() {
        auto _density_ratio = state.density_ratio;
        auto _viscosity_ratio = state.viscosity_ratio;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_xp = state.viscosity_ratio_xp;
        auto _viscosity_ratio_yp = state.viscosity_ratio_yp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _density_ratio_xp(domain_idx) = 1.0;
                    _density_ratio_yp(domain_idx) = 1.0;
                    _density_ratio_zp(domain_idx) = 1.0;
                    _viscosity_ratio_xp(domain_idx) = 1.0;
                    _viscosity_ratio_yp(domain_idx) = 1.0;
                    _viscosity_ratio_zp(domain_idx) = 1.0;
                    return;
                }
                
                // X-direction
                Ordinal active_right = _neighbor_right(i);
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    _density_ratio_xp(domain_idx) = 0.5 * (_density_ratio(domain_idx) + _density_ratio(idx_right));
                    _viscosity_ratio_xp(domain_idx) = 0.5 * (_viscosity_ratio(domain_idx) + _viscosity_ratio(idx_right));
                } else {
                    _density_ratio_xp(domain_idx) = _density_ratio(domain_idx);
                    _viscosity_ratio_xp(domain_idx) = _viscosity_ratio(domain_idx);
                }
                
                // Y-direction
                Ordinal active_front = _neighbor_front(i);
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    _density_ratio_yp(domain_idx) = 0.5 * (_density_ratio(domain_idx) + _density_ratio(idx_front));
                    _viscosity_ratio_yp(domain_idx) = 0.5 * (_viscosity_ratio(domain_idx) + _viscosity_ratio(idx_front));
                } else {
                    _density_ratio_yp(domain_idx) = _density_ratio(domain_idx);
                    _viscosity_ratio_yp(domain_idx) = _viscosity_ratio(domain_idx);
                }
                
                // Z-direction
                Ordinal active_top = _neighbor_top(i);
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    _density_ratio_zp(domain_idx) = 0.5 * (_density_ratio(domain_idx) + _density_ratio(idx_top));
                    _viscosity_ratio_zp(domain_idx) = 0.5 * (_viscosity_ratio(domain_idx) + _viscosity_ratio(idx_top));
                } else {
                    _density_ratio_zp(domain_idx) = _density_ratio(domain_idx);
                    _viscosity_ratio_zp(domain_idx) = _viscosity_ratio(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // ASSEMBLE LINEAR SYSTEM (Matrix Coefficients and RHS)
    // ========================================================================
    void assemble_linear_system() {
        compute_matrix_coefficients();
        compute_rhs();
    }
    
    // Compute matrix coefficients (Gxp, Gxm, Gyp, Gym, Gzp, Gzm, Gct)
    void compute_matrix_coefficients() {
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_zp = state.matrix_zp;
        auto _matrix_zm = state.matrix_zm;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_xp = state.viscosity_ratio_xp;
        auto _viscosity_ratio_yp = state.viscosity_ratio_yp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _water_content = state.water_content;
        auto _water_content_old = state.water_content_old;
        auto _water_content_sat = state.water_content_sat;
        auto _density_ratio = state.density_ratio;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        // Get layer thicknesses from domain
        auto _dz_layers = domain.dz_layers;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _matrix_diag(domain_idx) = 1.0;
                    _matrix_xp(domain_idx) = 0.0; _matrix_xm(domain_idx) = 0.0;
                    _matrix_yp(domain_idx) = 0.0; _matrix_ym(domain_idx) = 0.0;
                    _matrix_zp(domain_idx) = 0.0; _matrix_zm(domain_idx) = 0.0;
                    return;
                }
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal active_bottom = _neighbor_bottom(i);
                Ordinal active_top = _neighbor_top(i);
                
                // Get layer index (k)
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar dz_avg = dz; // Simplified: use current layer thickness
                
                // Compute specific moisture capacity
                Scalar ch = compute_ch(domain_idx);
                
                // Cell volume
                Scalar V = domain.dx * domain.dy * dz;
                
                // X-direction coefficients
                _matrix_xp(domain_idx) = 0.0;
                if (active_right >= 0) {
                    Scalar Kx = _conductivity_x(domain_idx);
                    Scalar r_rho = _density_ratio_xp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_xp(domain_idx);
                    Scalar Ax = domain.dy * dz; // Face area
                    _matrix_xp(domain_idx) = -Kx * dt * r_rho * r_visc * Ax / (domain.dx * domain.dx);
                }
                
                _matrix_xm(domain_idx) = 0.0;
                if (active_left >= 0) {
                    Ordinal idx_left = _active_to_domain(active_left);
                    Scalar Kx = _conductivity_x(idx_left);
                    Scalar r_rho = _density_ratio_xp(idx_left);
                    Scalar r_visc = _viscosity_ratio_xp(idx_left);
                    Scalar Ax = domain.dy * dz;
                    _matrix_xm(domain_idx) = -Kx * dt * r_rho * r_visc * Ax / (domain.dx * domain.dx);
                }
                
                // Y-direction coefficients
                _matrix_yp(domain_idx) = 0.0;
                if (active_front >= 0) {
                    Scalar Ky = _conductivity_y(domain_idx);
                    Scalar r_rho = _density_ratio_yp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_yp(domain_idx);
                    Scalar Ay = domain.dx * dz;
                    _matrix_yp(domain_idx) = -Ky * dt * r_rho * r_visc * Ay / (domain.dy * domain.dy);
                }
                
                _matrix_ym(domain_idx) = 0.0;
                if (active_back >= 0) {
                    Ordinal idx_back = _active_to_domain(active_back);
                    Scalar Ky = _conductivity_y(idx_back);
                    Scalar r_rho = _density_ratio_yp(idx_back);
                    Scalar r_visc = _viscosity_ratio_yp(idx_back);
                    Scalar Ay = domain.dx * dz;
                    _matrix_ym(domain_idx) = -Ky * dt * r_rho * r_visc * Ay / (domain.dy * domain.dy);
                }
                
                // Z-direction coefficients
                _matrix_zp(domain_idx) = 0.0;
                if (active_top >= 0) {
                    Scalar Kz = _conductivity_z(domain_idx);
                    Scalar r_rho = _density_ratio_zp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_zp(domain_idx);
                    Scalar Az = domain.dx * domain.dy;
                    Scalar dz_top = (k < domain.nz - 1) ? _dz_layers(k + 1) : dz;
                    Scalar dz_avg_z = 0.5 * (dz + dz_top);
                    _matrix_zp(domain_idx) = -Kz * dt * r_rho * r_visc * V / (dz * dz_avg_z);
                }
                
                _matrix_zm(domain_idx) = 0.0;
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    Scalar Kz = _conductivity_z(idx_bottom);
                    Scalar r_rho = _density_ratio_zp(idx_bottom);
                    Scalar r_visc = _viscosity_ratio_zp(idx_bottom);
                    Scalar Az = domain.dx * domain.dy;
                    Ordinal k_bottom = _coord_k(active_bottom);
                    Scalar dz_bottom = _dz_layers(k_bottom);
                    Scalar dz_avg_z = 0.5 * (dz + dz_bottom);
                    _matrix_zm(domain_idx) = -Kz * dt * r_rho * r_visc * V / (dz * dz_avg_z);
                }
                
                // Diagonal coefficient
                Scalar ch_term = (ch + Ss * _water_content_old(domain_idx) / _water_content_sat(domain_idx)) * 
                                 _density_ratio(domain_idx);
                _matrix_diag(domain_idx) = ch_term * V;
                _matrix_diag(domain_idx) -= (_matrix_xp(domain_idx) + _matrix_xm(domain_idx) + 
                                            _matrix_yp(domain_idx) + _matrix_ym(domain_idx) +
                                            _matrix_zp(domain_idx) + _matrix_zm(domain_idx));
            });
    }
    
    // Compute right-hand side vector
    void compute_rhs() {
        auto _matrix_rhs = state.matrix_rhs;
        auto _pressure = state.pressure;
        auto _pressure_old = state.pressure_old;
        auto _water_content = state.water_content;
        auto _water_content_old = state.water_content_old;
        auto _water_content_sat = state.water_content_sat;
        auto _density_ratio = state.density_ratio;
        auto _conductivity_z = state.conductivity_z;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        auto _dz_layers = domain.dz_layers;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _matrix_rhs(domain_idx) = _pressure_old(domain_idx);
                    return;
                }
                
                // Get layer index
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar V = domain.dx * domain.dy * dz;
                
                // Compute specific moisture capacity
                Scalar ch = compute_ch(domain_idx);
                
                // Base term: storage term
                Scalar ch_term = (ch + Ss * _water_content_old(domain_idx) / _water_content_sat(domain_idx)) * 
                                _density_ratio(domain_idx);
                _matrix_rhs(domain_idx) = ch_term * _pressure_old(domain_idx) * V;
                
                // Density-driven flow term (vertical)
                _matrix_rhs(domain_idx) -= dt * V * _conductivity_z(domain_idx) * 
                                          _density_ratio_zp(domain_idx) * _density_ratio_zp(domain_idx) * 
                                          _viscosity_ratio_zp(domain_idx) / dz;
                
                Ordinal active_bottom = _neighbor_bottom(i);
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    Ordinal k_bottom = _coord_k(active_bottom);
                    Scalar dz_bottom = _dz_layers(k_bottom);
                    _matrix_rhs(domain_idx) += dt * V * _conductivity_z(idx_bottom) * 
                                              _density_ratio_zp(idx_bottom) * _density_ratio_zp(idx_bottom) * 
                                              _viscosity_ratio_zp(idx_bottom) / dz_bottom;
                }
                
                // Boundary conditions and source/sink terms will be added here later
            });
    }
    
    // ========================================================================
    // SOLVE LINEAR SYSTEM
    // ========================================================================
    void solve_linear_system() {
        // Initialize solution vector with current head as initial guess
        auto _pressure = state.pressure;
        auto _active_to_domain = active_mesh.active_to_domain;
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                solution(i) = _pressure(domain_idx);
            });
        
        pcg_solver->solve_3d(active_mesh,
                            state.matrix_diag,
                            state.matrix_xp,
                            state.matrix_xm,
                            state.matrix_yp,
                            state.matrix_ym,
                            state.matrix_zp,
                            state.matrix_zm,
                            state.matrix_rhs,
                            solution);
        
        // Copy solution back to state.pressure (h_new)
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                _pressure(domain_idx) = solution(i);
            });
    }
    
    // ========================================================================
    // APPLY BOUNDARY CONDITIONS TO MATRIX
    // ========================================================================
    void apply_boundary_conditions_to_matrix(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_zp = state.matrix_zp;
        auto _matrix_zm = state.matrix_zm;
        auto _matrix_rhs = state.matrix_rhs;
        
        const auto& bcs = bc_manager_->get_boundary_conditions();
        
        // Create host mirrors for boundary condition application
        auto h_matrix_diag = Kokkos::create_mirror_view(_matrix_diag);
        auto h_matrix_xp = Kokkos::create_mirror_view(_matrix_xp);
        auto h_matrix_xm = Kokkos::create_mirror_view(_matrix_xm);
        auto h_matrix_yp = Kokkos::create_mirror_view(_matrix_yp);
        auto h_matrix_ym = Kokkos::create_mirror_view(_matrix_ym);
        auto h_matrix_zp = Kokkos::create_mirror_view(_matrix_zp);
        auto h_matrix_zm = Kokkos::create_mirror_view(_matrix_zm);
        auto h_matrix_rhs = Kokkos::create_mirror_view(_matrix_rhs);
        
        Kokkos::deep_copy(h_matrix_diag, _matrix_diag);
        Kokkos::deep_copy(h_matrix_xp, _matrix_xp);
        Kokkos::deep_copy(h_matrix_xm, _matrix_xm);
        Kokkos::deep_copy(h_matrix_yp, _matrix_yp);
        Kokkos::deep_copy(h_matrix_ym, _matrix_ym);
        Kokkos::deep_copy(h_matrix_zp, _matrix_zp);
        Kokkos::deep_copy(h_matrix_zm, _matrix_zm);
        Kokkos::deep_copy(h_matrix_rhs, _matrix_rhs);
        
        // Apply each boundary condition
        for (const auto& bc : bcs) {
            Scalar bc_value = bc.get_value(current_time);
            
            for (Ordinal cell_idx : bc.cell_indices) {
                if (bc.type == GwBcType::FIXED_HEAD) {
                    // Dirichlet BC: prescribed hydraulic head
                    h_matrix_diag(cell_idx) = 1.0;
                    h_matrix_xp(cell_idx) = 0.0;
                    h_matrix_xm(cell_idx) = 0.0;
                    h_matrix_yp(cell_idx) = 0.0;
                    h_matrix_ym(cell_idx) = 0.0;
                    h_matrix_zp(cell_idx) = 0.0;
                    h_matrix_zm(cell_idx) = 0.0;
                    h_matrix_rhs(cell_idx) = bc_value;
                    
                } else if (bc.type == GwBcType::FIXED_FLUX) {
                    // Neumann BC: prescribed flux (add to RHS)
                    Scalar flux_volume = bc_value * dt;
                    h_matrix_rhs(cell_idx) += flux_volume;
                }
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(_matrix_diag, h_matrix_diag);
        Kokkos::deep_copy(_matrix_xp, h_matrix_xp);
        Kokkos::deep_copy(_matrix_xm, h_matrix_xm);
        Kokkos::deep_copy(_matrix_yp, h_matrix_yp);
        Kokkos::deep_copy(_matrix_ym, h_matrix_ym);
        Kokkos::deep_copy(_matrix_zp, h_matrix_zp);
        Kokkos::deep_copy(_matrix_zm, h_matrix_zm);
        Kokkos::deep_copy(_matrix_rhs, h_matrix_rhs);
    }
    
    // ========================================================================
    // ENFORCE BOUNDARY CONDITIONS ON SOLUTION
    // ========================================================================
    void enforce_boundary_conditions_on_solution(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _pressure = state.pressure;
        const auto& bcs = bc_manager_->get_boundary_conditions();
        
        auto h_pressure = Kokkos::create_mirror_view(_pressure);
        Kokkos::deep_copy(h_pressure, _pressure);
        
        for (const auto& bc : bcs) {
            Scalar bc_value = bc.get_value(current_time);
            
            for (Ordinal cell_idx : bc.cell_indices) {
                if (bc.type == GwBcType::FIXED_HEAD) {
                    h_pressure(cell_idx) = bc_value;
                }
            }
        }
        
        Kokkos::deep_copy(_pressure, h_pressure);
    }
    
    // ========================================================================
    // COMPUTE GROUNDWATER FLUX
    // ========================================================================
    void compute_groundwater_flux() {
        // Compute Darcy fluxes in x, y, z directions
        // This is a simplified version - full implementation would use darcy_flux function
        auto _pressure = state.pressure;
        auto _conductivity_x = state.conductivity_x;
        auto _conductivity_y = state.conductivity_y;
        auto _conductivity_z = state.conductivity_z;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _viscosity_ratio_xp = state.viscosity_ratio_xp;
        auto _viscosity_ratio_yp = state.viscosity_ratio_yp;
        auto _viscosity_ratio_zp = state.viscosity_ratio_zp;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _active_mask_3d = domain.active_mask_3d;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_top = active_mesh.neighbor_top;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        auto _dz_layers = domain.dz_layers;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _flux_x(domain_idx) = 0.0;
                    _flux_y(domain_idx) = 0.0;
                    _flux_z(domain_idx) = 0.0;
                    return;
                }
                
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                
                // X-direction flux (Darcy's law)
                Ordinal active_right = _neighbor_right(i);
                if (active_right >= 0) {
                    Ordinal idx_right = _active_to_domain(active_right);
                    Scalar dh = _pressure(idx_right) - _pressure(domain_idx);
                    Scalar Kx = _conductivity_x(domain_idx);
                    Scalar r_rho = _density_ratio_xp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_xp(domain_idx);
                    Scalar Ax = domain.dy * dz;
                    _flux_x(domain_idx) = -Kx * r_visc * r_rho * Ax * dh / domain.dx;
                } else {
                    _flux_x(domain_idx) = 0.0;
                }
                
                // Y-direction flux
                Ordinal active_front = _neighbor_front(i);
                if (active_front >= 0) {
                    Ordinal idx_front = _active_to_domain(active_front);
                    Scalar dh = _pressure(idx_front) - _pressure(domain_idx);
                    Scalar Ky = _conductivity_y(domain_idx);
                    Scalar r_rho = _density_ratio_yp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_yp(domain_idx);
                    Scalar Ay = domain.dx * dz;
                    _flux_y(domain_idx) = -Ky * r_visc * r_rho * Ay * dh / domain.dy;
                } else {
                    _flux_y(domain_idx) = 0.0;
                }
                
                // Z-direction flux (includes density term)
                Ordinal active_top = _neighbor_top(i);
                if (active_top >= 0) {
                    Ordinal idx_top = _active_to_domain(active_top);
                    Scalar dh = _pressure(idx_top) - _pressure(domain_idx);
                    Scalar Kz = _conductivity_z(domain_idx);
                    Scalar r_rho = _density_ratio_zp(domain_idx);
                    Scalar r_visc = _viscosity_ratio_zp(domain_idx);
                    Scalar Az = domain.dx * domain.dy;
                    Scalar dz_top = (k < domain.nz - 1) ? _dz_layers(k + 1) : dz;
                    Scalar dz_avg = 0.5 * (dz + dz_top);
                    _flux_z(domain_idx) = -Kz * r_visc * r_rho * Az * (dh / dz_avg - r_rho);
                } else {
                    _flux_z(domain_idx) = 0.0;
                }
            });
    }
    
    // ========================================================================
    // CHECK AVAILABLE ROOM (Pore Space)
    // ========================================================================
    void check_room() {
        auto _water_content = state.water_content;
        auto _water_content_sat = state.water_content_sat;
        auto _room = state.room;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.dz_layers;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    _room(domain_idx) = 0.0;
                    return;
                }
                
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar V = domain.dx * domain.dy * dz;
                _room(domain_idx) = (_water_content_sat(domain_idx) - _water_content(domain_idx)) * V;
            });
    }
    
    // ========================================================================
    // UPDATE WATER CONTENT (From Fluxes)
    // ========================================================================
    void update_water_content() {
        auto _water_content = state.water_content;
        auto _water_content_old = state.water_content_old;
        auto _water_content_sat = state.water_content_sat;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _flux_z = state.flux_z;
        auto _density_ratio = state.density_ratio;
        auto _density_ratio_xp = state.density_ratio_xp;
        auto _density_ratio_yp = state.density_ratio_yp;
        auto _density_ratio_zp = state.density_ratio_zp;
        auto _pressure = state.pressure;
        auto _pressure_old = state.pressure_old;
        auto _active_mask_3d = domain.active_mask_3d;
        auto _dz_layers = domain.dz_layers;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_bottom = active_mesh.neighbor_bottom;
        auto _active_to_domain = active_mesh.active_to_domain;
        auto _coord_k = active_mesh.coord_k;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                if (_active_mask_3d(domain_idx) == 0) {
                    return;
                }
                
                Ordinal k = _coord_k(i);
                Scalar dz = _dz_layers(k);
                Scalar V = domain.dx * domain.dy * dz;
                
                // Compute specific moisture capacity
                Scalar ch = compute_ch(domain_idx);
                
                // Storage coefficient
                Scalar coeff = V * (_density_ratio(domain_idx) + 
                                   _density_ratio(domain_idx) * Ss * 
                                   (_pressure(domain_idx) - _pressure_old(domain_idx)) / 
                                   _water_content_sat(domain_idx));
                
                // Flux divergence
                Scalar dqx = 0.0;
                Ordinal active_left = _neighbor_left(i);
                if (active_left >= 0) {
                    Ordinal idx_left = _active_to_domain(active_left);
                    dqx = dt * (_flux_x(domain_idx) * _density_ratio_xp(domain_idx) - 
                               _flux_x(idx_left) * _density_ratio_xp(idx_left));
                } else {
                    dqx = dt * _flux_x(domain_idx) * _density_ratio_xp(domain_idx);
                }
                
                Scalar dqy = 0.0;
                Ordinal active_back = _neighbor_back(i);
                if (active_back >= 0) {
                    Ordinal idx_back = _active_to_domain(active_back);
                    dqy = dt * (_flux_y(domain_idx) * _density_ratio_yp(domain_idx) - 
                               _flux_y(idx_back) * _density_ratio_yp(idx_back));
                } else {
                    dqy = dt * _flux_y(domain_idx) * _density_ratio_yp(domain_idx);
                }
                
                Scalar dqz = 0.0;
                Ordinal active_bottom = _neighbor_bottom(i);
                if (active_bottom >= 0) {
                    Ordinal idx_bottom = _active_to_domain(active_bottom);
                    dqz = dt * (_flux_z(domain_idx) * _density_ratio_zp(domain_idx) - 
                               _flux_z(idx_bottom) * _density_ratio_zp(idx_bottom));
                } else {
                    dqz = dt * _flux_z(domain_idx) * _density_ratio_zp(domain_idx);
                }
                
                // Update water content
                _water_content(domain_idx) = (_water_content_old(domain_idx) * _density_ratio(domain_idx) * V + 
                                             dqx + dqy + dqz) / coeff;
                
                // Clamp to valid range
                if (_water_content(domain_idx) > _water_content_sat(domain_idx)) {
                    _water_content(domain_idx) = _water_content_sat(domain_idx);
                }
                if (_water_content(domain_idx) < state.water_content_res(domain_idx)) {
                    _water_content(domain_idx) = state.water_content_res(domain_idx);
                }
            });
    }
};

#endif // FREHG_GROUNDWATER_SOLVER_HPP

