#ifndef FREHG_SHALLOW_WATER_SOLVER_HPP
#define FREHG_SHALLOW_WATER_SOLVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "PCGSolver.hpp"
#include "BoundaryConditions.hpp"
#include "SourceSinkTerms.hpp"
#include <memory>
#include <cmath>

// ============================================================================
//                      SHALLOW WATER EQUATION SOLVER
// ============================================================================
// Solves the 2D depth-integrated Navier-Stokes equations (shallow water)
// Uses Kokkos for parallelization and PCG for linear system solution

class ShallowWaterSolver {
public:
    // --- Domain and Mesh ---
    SwDomain& domain;
    ActiveCellMesh& active_mesh;
    
    // --- State Variables ---
    SwStateVariables& state;
    
    // --- Physical Parameters ---
    Scalar dt;                  // Time step
    Scalar grav;                // Gravitational acceleration
    Scalar manning_n;           // Manning's roughness coefficient
    Scalar visc_x;              // Horizontal eddy viscosity (x-direction)
    Scalar visc_y;              // Horizontal eddy viscosity (y-direction)
    Scalar min_depth;           // Minimum depth threshold
    Scalar water_threshold;     // Water depth threshold for wetting (wtfh)
    
    // --- Linear Solver ---
    std::unique_ptr<PCGSolver2D> pcg_solver;
    PreconditionerType precond_type;
    Scalar solver_tolerance;
    Ordinal max_solver_iterations;
    
    // --- Solution Vector (active cells only) ---
    View1D<Scalar> solution;    // Solution vector for PCG (eta_new in active space)
    
    // --- Boundary Conditions ---
    SwBoundaryConditionManager* bc_manager_;  // Pointer to boundary condition manager
    
    // --- Source/Sink Terms ---
    SwSourceSinkManager* source_sink_manager_;  // Pointer to source/sink manager
    
    // Constructor
    ShallowWaterSolver(SwDomain& _domain,
                      ActiveCellMesh& _active_mesh,
                      SwStateVariables& _state,
                      Scalar _dt,
                      Scalar _grav = Constants::g,
                      Scalar _manning_n = 0.03,
                      Scalar _visc_x = 0.0,
                      Scalar _visc_y = 0.0,
                      Scalar _min_depth = 1.0e-6,
                      Scalar _water_threshold = 0.01,
                      PreconditionerType _precond = PreconditionerType::JACOBI,
                      Scalar _solver_tol = 1.0e-8,
                      Ordinal _max_iter = 10000,
                      SwBoundaryConditionManager* _bc_manager = nullptr,
                      SwSourceSinkManager* _source_sink_manager = nullptr)
        : domain(_domain),
          active_mesh(_active_mesh),
          state(_state),
          dt(_dt),
          grav(_grav),
          manning_n(_manning_n),
          visc_x(_visc_x),
          visc_y(_visc_y),
          min_depth(_min_depth),
          water_threshold(_water_threshold),
          precond_type(_precond),
          solver_tolerance(_solver_tol),
          max_solver_iterations(_max_iter),
          bc_manager_(_bc_manager),
          source_sink_manager_(_source_sink_manager) {
        
        // Create PCG solver
        pcg_solver = std::make_unique<PCGSolver2D>(
            active_mesh.num_active, precond_type, solver_tolerance, max_solver_iterations);
        
        // Allocate solution vector (active cells only)
        solution = View1D<Scalar>("solution", active_mesh.num_active);
    }
    
    // Set boundary condition manager
    void set_boundary_condition_manager(SwBoundaryConditionManager* bc_manager) {
        bc_manager_ = bc_manager;
    }
    
    // Set source/sink manager
    void set_source_sink_manager(SwSourceSinkManager* source_sink_manager) {
        source_sink_manager_ = source_sink_manager;
    }
    
    // ========================================================================
    // MAIN SOLVER: Solve for surface elevation (eta)
    // ========================================================================
    void solve(Scalar current_time = 0.0) {
        // Step 1: Save old values
        state.update_old_values();
        
        // Step 2: Compute momentum source terms (Ex, Ey, Dx, Dy)
        compute_momentum_source();
        
        // Step 3: Compute matrix RHS
        compute_rhs();
        
        // Step 4: Compute matrix coefficients
        compute_matrix_coefficients();
        
        // Step 5: Apply boundary conditions to matrix
        apply_boundary_conditions_to_matrix(current_time);
        
        // Step 6: Solve linear system
        solve_linear_system();
        
        // Step 7: Map solution back to domain
        map_solution_to_domain();
        
        // Step 8: Enforce boundary conditions on solution
        enforce_boundary_conditions_on_solution(current_time);
        
        // Step 9: CFL limiter
        apply_cfl_limiter();
        
        // Step 10: Apply source/sink terms
        apply_source_sink_terms(current_time);
        
        // Step 11: Update depth
        update_depth();
    }
    
    // ========================================================================
    // VELOCITY UPDATE: Update velocities from pressure gradient
    // ========================================================================
    void update_velocity() {
        // Step 1: Update depth (if needed)
        update_depth();
        
        // Step 2: Update geometry (volumes and areas)
        update_geometry();
        
        // Step 3: Update drag coefficients
        update_drag_coefficients();
        
        // Step 4: Update velocity from pressure gradient
        compute_velocity_from_pressure();
        
        // Step 5: Apply velocity limiters
        apply_velocity_limiters();
        
        // Step 6: Update fluxes and CFL numbers
        update_fluxes_and_cfl();
        
        // Step 7: Interpolate velocities to cell centers
        interpolate_velocities();
    }
    
private:
    // ========================================================================
    // MOMENTUM SOURCE TERMS
    // ========================================================================
    void compute_momentum_source() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _velocity_x_at_y = state.velocity_x_at_y;
        auto _velocity_y_at_x = state.velocity_y_at_x;
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        auto _drag_coef_x = state.drag_coef_x;
        auto _drag_coef_y = state.drag_coef_y;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _area_top_x = state.area_top_x;
        auto _area_top_y = state.area_top_y;
        auto _cfl_x = state.cfl_x;
        auto _cfl_y = state.cfl_y;
        auto _depth = state.depth;
        auto _depth_x = state.depth_x;
        auto _depth_y = state.depth_y;
        auto _pressure = state.pressure;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                // Get domain index
                Ordinal domain_idx = _active_to_domain(i);
                
                // Get neighbor active indices
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_front = _neighbor_front(i);
                
                // Get neighbor domain indices (neighbors are stored as active indices)
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                
                // Advection terms (upwind scheme)
                Scalar adv_x = 0.0;
                Scalar adv_y = 0.0;
                
                // X-direction advection
                if (active_left >= 0 && active_right >= 0) {
                    Scalar u_here = _velocity_x(domain_idx);
                    Scalar u_left = (i_left >= 0) ? _velocity_x(i_left) : u_here;
                    Scalar u_right = (i_right >= 0) ? _velocity_x(i_right) : u_here;
                    Scalar v_interp = _velocity_y_at_x(domain_idx);
                    Scalar u_back = (i_back >= 0) ? _velocity_x(i_back) : u_here;
                    Scalar u_front = (i_front >= 0) ? _velocity_x(i_front) : u_here;
                    
                    adv_x = (0.5 / domain.dx) * 
                        ((u_here + std::abs(u_here)) * (u_here - u_left) +
                         (u_here - std::abs(u_here)) * (u_right - u_here)) +
                        (0.5 / domain.dy) *
                        ((v_interp + std::abs(v_interp)) * (u_here - u_back) +
                         (v_interp - std::abs(v_interp)) * (u_front - u_here));
                }
                
                // Y-direction advection
                if (active_back >= 0 && active_front >= 0) {
                    Scalar v_here = _velocity_y(domain_idx);
                    Scalar v_back = (i_back >= 0) ? _velocity_y(i_back) : v_here;
                    Scalar v_front = (i_front >= 0) ? _velocity_y(i_front) : v_here;
                    Scalar u_interp = _velocity_x_at_y(domain_idx);
                    Scalar v_left = (i_left >= 0) ? _velocity_y(i_left) : v_here;
                    Scalar v_right = (i_right >= 0) ? _velocity_y(i_right) : v_here;
                    
                    adv_y = (0.5 / domain.dx) *
                        ((u_interp + std::abs(u_interp)) * (v_here - v_left) +
                         (u_interp - std::abs(u_interp)) * (v_right - v_here)) +
                        (0.5 / domain.dy) *
                        ((v_here + std::abs(v_here)) * (v_here - v_back) +
                         (v_here - std::abs(v_here)) * (v_front - v_here));
                }
                
                // CFL-based advection limiter
                if (_velocity_x(domain_idx) == 0.0 || _cfl_x(domain_idx) > 0.7) {
                    adv_x = 0.0;
                } else if (_cfl_x(domain_idx) > 0.5) {
                    adv_x *= (0.7 - _cfl_x(domain_idx)) / (0.7 - 0.5);
                }
                
                if (_velocity_y(domain_idx) == 0.0 || _cfl_y(domain_idx) > 0.7) {
                    adv_y = 0.0;
                } else if (_cfl_y(domain_idx) > 0.5) {
                    adv_y *= (0.7 - _cfl_y(domain_idx)) / (0.7 - 0.5);
                }
                
                // Diffusion terms
                Scalar diff_x = 0.0;
                Scalar diff_y = 0.0;
                
                if (_volume_x(domain_idx) > 0.0 && visc_x > 0.0) {
                    if (active_left >= 0 && active_right >= 0) {
                        diff_x = (visc_x / _volume_x(domain_idx) / domain.dx) *
                            (_area_x(domain_idx) * (_velocity_x(i_right) - _velocity_x(domain_idx)) -
                             _area_x(domain_idx) * (_velocity_x(domain_idx) - _velocity_x(i_left)));
                    }
                    if (active_back >= 0 && active_front >= 0) {
                        Scalar area_y_back = (i_back >= 0) ? _area_y(i_back) : 0.0;
                        diff_x += (visc_y / _volume_x(domain_idx) / domain.dy) *
                            (_area_y(domain_idx) * (_velocity_x(i_front) - _velocity_x(domain_idx)) -
                             area_y_back * (_velocity_x(domain_idx) - _velocity_x(i_back)));
                    }
                }
                
                if (_volume_y(domain_idx) > 0.0 && visc_y > 0.0) {
                    if (active_left >= 0 && active_right >= 0) {
                        Scalar area_x_left = (i_left >= 0) ? _area_x(i_left) : 0.0;
                        diff_y = (visc_x / _volume_y(domain_idx) / domain.dx) *
                            (_area_x(domain_idx) * (_velocity_y(i_right) - _velocity_y(domain_idx)) -
                             area_x_left * (_velocity_y(domain_idx) - _velocity_y(i_left)));
                    }
                    if (active_back >= 0 && active_front >= 0) {
                        diff_y += (visc_y / _volume_y(domain_idx) / domain.dy) *
                            (_area_y(domain_idx) * (_velocity_y(i_front) - _velocity_y(domain_idx)) -
                             _area_y(domain_idx) * (_velocity_y(domain_idx) - _velocity_y(i_back)));
                    }
                }
                
                // Drag terms
                Scalar fac_dx = 0.0;
                Scalar fac_dy = 0.0;
                Scalar vel_mag_x = std::sqrt(_velocity_x(domain_idx) * _velocity_x(domain_idx) + 
                                            _velocity_y_at_x(domain_idx) * _velocity_y_at_x(domain_idx));
                Scalar vel_mag_y = std::sqrt(_velocity_x_at_y(domain_idx) * _velocity_x_at_y(domain_idx) + 
                                            _velocity_y(domain_idx) * _velocity_y(domain_idx));
                
                if (_volume_x(domain_idx) > 0.0) {
                    fac_dx = _area_top_x(domain_idx) / _volume_x(domain_idx);
                }
                if (_volume_y(domain_idx) > 0.0) {
                    fac_dy = _area_top_y(domain_idx) / _volume_y(domain_idx);
                }
                
                // Standard momentum source
                _drag_factor_x(domain_idx) = 1.0 / (0.5 * dt * _drag_coef_x(domain_idx) * vel_mag_x * fac_dx + 1.0);
                _drag_factor_y(domain_idx) = 1.0 / (0.5 * dt * _drag_coef_y(domain_idx) * vel_mag_y * fac_dy + 1.0);
                
                _momentum_x(domain_idx) = _velocity_x(domain_idx) + dt * (diff_x - adv_x);
                _momentum_y(domain_idx) = _velocity_y(domain_idx) + dt * (diff_y - adv_y);
                
                // Apply drag factor
                _momentum_x(domain_idx) *= _drag_factor_x(domain_idx);
                _momentum_y(domain_idx) *= _drag_factor_y(domain_idx);
            });
    }
    
    // ========================================================================
    // COMPUTE RHS
    // ========================================================================
    void compute_rhs() {
        auto _pressure = state.pressure;
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _area_top = state.area_top;
        auto _matrix_rhs = state.matrix_rhs;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                // Get neighbors
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                
                // RHS = eta * Asz - dt * div(Ex, Ey)
                Scalar div_flux = _area_x(domain_idx) * _momentum_x(domain_idx);
                if (i_left >= 0) {
                    div_flux -= _area_x(i_left) * _momentum_x(i_left);
                }
                
                div_flux += _area_y(domain_idx) * _momentum_y(domain_idx);
                if (i_back >= 0) {
                    div_flux -= _area_y(i_back) * _momentum_y(i_back);
                }
                
                _matrix_rhs(domain_idx) = _pressure(domain_idx) * _area_top(domain_idx) - dt * div_flux;
                
                // Inflow source terms would be added here (assumed available)
                // if (has_inflow(domain_idx)) {
                //     _matrix_rhs(domain_idx) += inflow_rate * dt;
                // }
            });
    }
    
    // ========================================================================
    // COMPUTE MATRIX COEFFICIENTS
    // ========================================================================
    void compute_matrix_coefficients() {
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_rhs = state.matrix_rhs;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _area_top = state.area_top;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        auto _depth = state.depth;
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _pressure = state.pressure;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        
        Scalar coef = grav * dt * dt;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                
                // Initialize coefficients
                _matrix_xp(domain_idx) = 0.0;
                _matrix_xm(domain_idx) = 0.0;
                _matrix_yp(domain_idx) = 0.0;
                _matrix_ym(domain_idx) = 0.0;
                
                // X+ coefficient
                if (_volume_x(domain_idx) > 0.0) {
                    _matrix_xp(domain_idx) = coef * _area_x(domain_idx) * _area_x(domain_idx) * 
                                            _drag_factor_x(domain_idx) / _volume_x(domain_idx);
                }
                
                // X- coefficient (from left neighbor)
                if (i_left >= 0 && _volume_x(i_left) > 0.0) {
                    _matrix_xm(domain_idx) = coef * _area_x(i_left) * _area_x(i_left) * 
                                            _drag_factor_x(i_left) / _volume_x(i_left);
                }
                
                // Y+ coefficient
                if (_volume_y(domain_idx) > 0.0) {
                    _matrix_yp(domain_idx) = coef * _area_y(domain_idx) * _area_y(domain_idx) * 
                                            _drag_factor_y(domain_idx) / _volume_y(domain_idx);
                }
                
                // Y- coefficient (from back neighbor)
                if (i_back >= 0 && _volume_y(i_back) > 0.0) {
                    _matrix_ym(domain_idx) = coef * _area_y(i_back) * _area_y(i_back) * 
                                            _drag_factor_y(i_back) / _volume_y(i_back);
                }
                
                // Diagonal
                _matrix_diag(domain_idx) = _area_top(domain_idx) + 
                                          _matrix_xp(domain_idx) + _matrix_xm(domain_idx) +
                                          _matrix_yp(domain_idx) + _matrix_ym(domain_idx);
                
                // Avoid singularity
                if (_depth(domain_idx) == 0.0) {
                    _matrix_diag(domain_idx) = domain.dx * domain.dy;
                    _matrix_rhs(domain_idx) = _pressure(domain_idx) * domain.dx * domain.dy;
                    
                    if (_velocity_x(domain_idx) == 0.0 && 
                        (i_left < 0 || _velocity_x(i_left) == 0.0) &&
                        _velocity_y(domain_idx) == 0.0 &&
                        (i_back < 0 || _velocity_y(i_back) == 0.0)) {
                        _matrix_xp(domain_idx) = 0.0;
                        _matrix_xm(domain_idx) = 0.0;
                        _matrix_yp(domain_idx) = 0.0;
                        _matrix_ym(domain_idx) = 0.0;
                    }
                } else if (_matrix_diag(domain_idx) == 0.0) {
                    _matrix_diag(domain_idx) = domain.dx * domain.dy;
                    _matrix_rhs(domain_idx) = _pressure(domain_idx) * domain.dx * domain.dy;
                }
                
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
        auto _matrix_rhs = state.matrix_rhs;
        auto _bottom = state.bottom;
        
        const auto& bcs = bc_manager_->get_boundary_conditions();
        
        // Create host mirrors for boundary condition application
        auto h_matrix_diag = Kokkos::create_mirror_view(_matrix_diag);
        auto h_matrix_xp = Kokkos::create_mirror_view(_matrix_xp);
        auto h_matrix_xm = Kokkos::create_mirror_view(_matrix_xm);
        auto h_matrix_yp = Kokkos::create_mirror_view(_matrix_yp);
        auto h_matrix_ym = Kokkos::create_mirror_view(_matrix_ym);
        auto h_matrix_rhs = Kokkos::create_mirror_view(_matrix_rhs);
        auto h_bottom = Kokkos::create_mirror_view(_bottom);
        
        Kokkos::deep_copy(h_matrix_diag, _matrix_diag);
        Kokkos::deep_copy(h_matrix_xp, _matrix_xp);
        Kokkos::deep_copy(h_matrix_xm, _matrix_xm);
        Kokkos::deep_copy(h_matrix_yp, _matrix_yp);
        Kokkos::deep_copy(h_matrix_ym, _matrix_ym);
        Kokkos::deep_copy(h_matrix_rhs, _matrix_rhs);
        Kokkos::deep_copy(h_bottom, _bottom);
        
        // Apply each boundary condition
        for (const auto& bc : bcs) {
            Scalar bc_value = bc.get_value(current_time);
            
            for (Ordinal cell_idx : bc.cell_indices) {
                if (bc.type == SwBcType::FREE_SURFACE_ELEVATION) {
                    // Dirichlet BC: prescribed surface elevation
                    // Set diagonal = 1.0, off-diagonals = 0.0, RHS = BC value
                    h_matrix_diag(cell_idx) = 1.0;
                    h_matrix_xp(cell_idx) = 0.0;
                    h_matrix_xm(cell_idx) = 0.0;
                    h_matrix_yp(cell_idx) = 0.0;
                    h_matrix_ym(cell_idx) = 0.0;
                    h_matrix_rhs(cell_idx) = bc_value;
                    
                } else if (bc.type == SwBcType::WATER_DEPTH) {
                    // Dirichlet BC: prescribed water depth
                    // Convert depth to surface elevation: eta = bottom + depth
                    Scalar eta_value = h_bottom(cell_idx) + bc_value;
                    h_matrix_diag(cell_idx) = 1.0;
                    h_matrix_xp(cell_idx) = 0.0;
                    h_matrix_xm(cell_idx) = 0.0;
                    h_matrix_yp(cell_idx) = 0.0;
                    h_matrix_ym(cell_idx) = 0.0;
                    h_matrix_rhs(cell_idx) = eta_value;
                    
                } else if (bc.type == SwBcType::FLOW_RATE) {
                    // Flow rate BC: add to RHS as source term
                    // Flow rate is discharge per unit width, convert to volume flux
                    // For now, treat as a source term in the RHS
                    // This is a simplified treatment - full implementation would
                    // modify flux terms at boundaries
                    Scalar source_flux = bc_value * dt;  // Volume added per time step
                    h_matrix_rhs(cell_idx) += source_flux;
                    
                } else if (bc.type == SwBcType::FREE_OUTFLOW) {
                    // Free outflow: zero gradient (Neumann BC)
                    // Remove boundary connection by zeroing off-diagonal coefficients
                    // This allows free flow out of the domain
                    // The diagonal remains unchanged, allowing natural outflow
                    // For cells at domain edge, we zero the appropriate off-diagonal
                    // For now, we'll keep the matrix as-is and let natural outflow occur
                    // This could be refined to explicitly set gradient = 0
                }
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(_matrix_diag, h_matrix_diag);
        Kokkos::deep_copy(_matrix_xp, h_matrix_xp);
        Kokkos::deep_copy(_matrix_xm, h_matrix_xm);
        Kokkos::deep_copy(_matrix_yp, h_matrix_yp);
        Kokkos::deep_copy(_matrix_ym, h_matrix_ym);
        Kokkos::deep_copy(_matrix_rhs, h_matrix_rhs);
    }
    
    // ========================================================================
    // ENFORCE BOUNDARY CONDITIONS ON SOLUTION
    // ========================================================================
    void enforce_boundary_conditions_on_solution(Scalar current_time) {
        if (!bc_manager_) return;
        
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        
        const auto& bcs = bc_manager_->get_boundary_conditions();
        
        // Create host mirrors
        auto h_pressure = Kokkos::create_mirror_view(_pressure);
        auto h_bottom = Kokkos::create_mirror_view(_bottom);
        
        Kokkos::deep_copy(h_pressure, _pressure);
        Kokkos::deep_copy(h_bottom, _bottom);
        
        // Apply each boundary condition
        for (const auto& bc : bcs) {
            Scalar bc_value = bc.get_value(current_time);
            
            for (Ordinal cell_idx : bc.cell_indices) {
                if (bc.type == SwBcType::FREE_SURFACE_ELEVATION) {
                    // Enforce prescribed surface elevation
                    h_pressure(cell_idx) = bc_value;
                    
                } else if (bc.type == SwBcType::WATER_DEPTH) {
                    // Enforce prescribed water depth
                    h_pressure(cell_idx) = h_bottom(cell_idx) + bc_value;
                    
                } else if (bc.type == SwBcType::FLOW_RATE) {
                    // Flow rate BC is handled in matrix, but ensure solution is valid
                    // (no additional enforcement needed here)
                    
                } else if (bc.type == SwBcType::FREE_OUTFLOW) {
                    // Free outflow: ensure solution is valid (no negative depth)
                    if (h_pressure(cell_idx) < h_bottom(cell_idx)) {
                        h_pressure(cell_idx) = h_bottom(cell_idx);
                    }
                }
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(_pressure, h_pressure);
    }
    
    // ========================================================================
    // APPLY SOURCE/SINK TERMS
    // ========================================================================
    void apply_source_sink_terms(Scalar current_time) {
        if (!source_sink_manager_) return;
        
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        auto _area_top = state.area_top;
        
        const auto& ss_terms = source_sink_manager_->get_source_sink_terms();
        
        // Create host mirrors
        auto h_pressure = Kokkos::create_mirror_view(_pressure);
        auto h_bottom = Kokkos::create_mirror_view(_bottom);
        auto h_area_top = Kokkos::create_mirror_view(_area_top);
        
        Kokkos::deep_copy(h_pressure, _pressure);
        Kokkos::deep_copy(h_bottom, _bottom);
        Kokkos::deep_copy(h_area_top, _area_top);
        
        // Apply each source/sink term
        for (const auto& ss : ss_terms) {
            Scalar ss_value = ss.get_value(current_time);
            
            for (Ordinal cell_idx : ss.cell_indices) {
                Scalar cell_area = h_area_top(cell_idx);
                if (cell_area <= 0.0) cell_area = domain.dx * domain.dy;
                
                if (ss.type == SourceSinkType::VOLUME_FLUX) {
                    // Volume flux (m³/s): add/remove volume, convert to depth change
                    Scalar volume_change = ss_value * dt;
                    Scalar depth_change = volume_change / cell_area;
                    h_pressure(cell_idx) += depth_change;
                    
                } else if (ss.type == SourceSinkType::DEPTH_RATE) {
                    // Depth rate (m/s): directly add/remove depth
                    Scalar depth_change = ss_value * dt;
                    h_pressure(cell_idx) += depth_change;
                    
                } else if (ss.type == SourceSinkType::MASS_FLUX) {
                    // Mass flux: for scalar transport only, not applied here
                    // This will be handled in scalar transport solver
                }
                
                // Ensure eta >= bottom
                if (h_pressure(cell_idx) < h_bottom(cell_idx)) {
                    h_pressure(cell_idx) = h_bottom(cell_idx);
                }
            }
        }
        
        // Copy back to device
        Kokkos::deep_copy(_pressure, h_pressure);
    }
    
    // ========================================================================
    // SOLVE LINEAR SYSTEM
    // ========================================================================
    void solve_linear_system() {
        // Create views for active cells only
        View1D<Scalar> diag_active("diag_active", active_mesh.num_active);
        View1D<Scalar> xp_active("xp_active", active_mesh.num_active);
        View1D<Scalar> xm_active("xm_active", active_mesh.num_active);
        View1D<Scalar> yp_active("yp_active", active_mesh.num_active);
        View1D<Scalar> ym_active("ym_active", active_mesh.num_active);
        View1D<Scalar> rhs_active("rhs_active", active_mesh.num_active);
        
        // Map domain variables to active cell space
        auto _matrix_diag = state.matrix_diag;
        auto _matrix_xp = state.matrix_xp;
        auto _matrix_xm = state.matrix_xm;
        auto _matrix_yp = state.matrix_yp;
        auto _matrix_ym = state.matrix_ym;
        auto _matrix_rhs = state.matrix_rhs;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                diag_active(i) = _matrix_diag(domain_idx);
                xp_active(i) = _matrix_xp(domain_idx);
                xm_active(i) = _matrix_xm(domain_idx);
                yp_active(i) = _matrix_yp(domain_idx);
                ym_active(i) = _matrix_ym(domain_idx);
                rhs_active(i) = _matrix_rhs(domain_idx);
            });
        
        // Solve using PCG
        pcg_solver->solve_2d(active_mesh,
                            diag_active, xp_active, xm_active, yp_active, ym_active,
                            rhs_active, solution);
    }
    
    // ========================================================================
    // MAP SOLUTION BACK TO DOMAIN
    // ========================================================================
    void map_solution_to_domain() {
        auto _pressure = state.pressure;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                _pressure(domain_idx) = solution(i);
            });
    }
    
    // ========================================================================
    // UPDATE DEPTH
    // ========================================================================
    void update_depth() {
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        auto _depth = state.depth;
        auto _depth_x = state.depth_x;
        auto _depth_y = state.depth_y;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Scalar h = _pressure(i) - _bottom(i);
                _depth(i) = (h > 0.0) ? h : 0.0;
                
                // Update face depths (simplified - would use neighbor info in full implementation)
                _depth_x(i) = _depth(i);
                _depth_y(i) = _depth(i);
            });
    }
    
    // ========================================================================
    // UPDATE GEOMETRY (Volumes and Areas)
    // ========================================================================
    void update_geometry() {
        auto _depth = state.depth;
        auto _area_top = state.area_top;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _volume = state.volume;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                // Top area
                if (_depth(i) > 0.0) {
                    _area_top(i) = domain.dx * domain.dy;
                } else {
                    _area_top(i) = 0.0;
                }
                
                // Face areas
                _area_x(i) = _depth(i) * domain.dy;
                _area_y(i) = _depth(i) * domain.dx;
                
                // Volumes
                _volume(i) = _depth(i) * domain.dx * domain.dy;
                
                // Staggered volumes (simplified - would average with neighbors)
                _volume_x(i) = _volume(i);
                _volume_y(i) = _volume(i);
            });
    }
    
    // ========================================================================
    // UPDATE DRAG COEFFICIENTS
    // ========================================================================
    void update_drag_coefficients() {
        auto _depth = state.depth;
        auto _volume = state.volume;
        auto _drag_coef_x = state.drag_coef_x;
        auto _drag_coef_y = state.drag_coef_y;
        
        Scalar coef = grav * manning_n * manning_n;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                if (_volume(i) > 0.0) {
                    Scalar eff_depth = _volume(i) / (domain.dx * domain.dy);
                    Scalar expo = (eff_depth < 0.1) ? (2.0 / 3.0) : (1.0 / 3.0);  // Thin layer model
                    Scalar cd = coef / std::pow(eff_depth, expo);
                    _drag_coef_x(i) = cd;
                    _drag_coef_y(i) = cd;
                } else {
                    _drag_coef_x(i) = 0.0;
                    _drag_coef_y(i) = 0.0;
                }
            });
    }
    
    // ========================================================================
    // COMPUTE VELOCITY FROM PRESSURE GRADIENT
    // ========================================================================
    void compute_velocity_from_pressure() {
        auto _pressure = state.pressure;
        auto _momentum_x = state.momentum_x;
        auto _momentum_y = state.momentum_y;
        auto _drag_factor_x = state.drag_factor_x;
        auto _drag_factor_y = state.drag_factor_y;
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _volume_x = state.volume_x;
        auto _volume_y = state.volume_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _depth = state.depth;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        
        Scalar coef = grav * dt;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                
                Scalar eff_hx = 0.0;
                Scalar eff_hy = 0.0;
                
                if (_volume_x(domain_idx) > 0.0) {
                    eff_hx = _area_x(domain_idx) / _volume_x(domain_idx);
                }
                if (_volume_y(domain_idx) > 0.0) {
                    eff_hy = _area_y(domain_idx) / _volume_y(domain_idx);
                }
                
                // Standard velocity update
                if (i_right >= 0) {
                    _velocity_x(domain_idx) = (_momentum_x(domain_idx) - 
                                               coef * eff_hx * (_pressure(i_right) - _pressure(domain_idx))) *
                                              _drag_factor_x(domain_idx);
                } else {
                    _velocity_x(domain_idx) = _momentum_x(domain_idx) * _drag_factor_x(domain_idx);
                }
                
                if (i_front >= 0) {
                    _velocity_y(domain_idx) = (_momentum_y(domain_idx) - 
                                              coef * eff_hy * (_pressure(i_front) - _pressure(domain_idx))) *
                                             _drag_factor_y(domain_idx);
                } else {
                    _velocity_y(domain_idx) = _momentum_y(domain_idx) * _drag_factor_y(domain_idx);
                }
            });
    }
    
    // ========================================================================
    // APPLY VELOCITY LIMITERS
    // ========================================================================
    void apply_velocity_limiters() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _depth = state.depth;
        auto _cfl_active = state.cfl_active;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                // Zero velocity when face area is too small
                if (_area_x(domain_idx) < water_threshold * domain.dy) {
                    _velocity_x(domain_idx) = 0.0;
                }
                if (_area_y(domain_idx) < water_threshold * domain.dx) {
                    _velocity_y(domain_idx) = 0.0;
                }
                
                // Zero velocity out of dry cells
                if (_depth(domain_idx) < water_threshold) {
                    if (_velocity_x(domain_idx) > 0.0) {
                        _velocity_x(domain_idx) = 0.0;
                    }
                    Ordinal active_left = _neighbor_left(i);
                    if (active_left >= 0) {
                        Ordinal i_left = _active_to_domain(active_left);
                        if (_velocity_x(i_left) < 0.0) {
                            _velocity_x(i_left) = 0.0;
                        }
                    }
                    
                    if (_velocity_y(domain_idx) > 0.0) {
                        _velocity_y(domain_idx) = 0.0;
                    }
                    Ordinal active_back = _neighbor_back(i);
                    if (active_back >= 0) {
                        Ordinal i_back = _active_to_domain(active_back);
                        if (_velocity_y(i_back) < 0.0) {
                            _velocity_y(i_back) = 0.0;
                        }
                    }
                }
                
                // CFL limiter
                if (_cfl_active(domain_idx) == 1) {
                    _velocity_x(domain_idx) = 0.0;
                    Ordinal active_left = _neighbor_left(i);
                    if (active_left >= 0) {
                        Ordinal i_left = _active_to_domain(active_left);
                        _velocity_x(i_left) = 0.0;
                    }
                    _velocity_y(domain_idx) = 0.0;
                    Ordinal active_back = _neighbor_back(i);
                    if (active_back >= 0) {
                        Ordinal i_back = _active_to_domain(active_back);
                        _velocity_y(i_back) = 0.0;
                    }
                    _cfl_active(domain_idx) = 0;
                }
            });
    }
    
    // ========================================================================
    // UPDATE FLUXES AND CFL NUMBERS
    // ========================================================================
    void update_fluxes_and_cfl() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _area_x = state.area_x;
        auto _area_y = state.area_y;
        auto _flux_x = state.flux_x;
        auto _flux_y = state.flux_y;
        auto _cfl_x = state.cfl_x;
        auto _cfl_y = state.cfl_y;
        
        Kokkos::parallel_for(RangePolicy(0, domain.num_cells_total),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                _flux_x(i) = _velocity_x(i) * _area_x(i);
                _flux_y(i) = _velocity_y(i) * _area_y(i);
                _cfl_x(i) = std::abs(_velocity_x(i) * dt / domain.dx);
                _cfl_y(i) = std::abs(_velocity_y(i) * dt / domain.dy);
            });
    }
    
    // ========================================================================
    // INTERPOLATE VELOCITIES
    // ========================================================================
    void interpolate_velocities() {
        auto _velocity_x = state.velocity_x;
        auto _velocity_y = state.velocity_y;
        auto _velocity_x_at_y = state.velocity_x_at_y;
        auto _velocity_y_at_x = state.velocity_y_at_x;
        
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Ordinal active_left = _neighbor_left(i);
                Ordinal active_back = _neighbor_back(i);
                Ordinal active_right = _neighbor_right(i);
                Ordinal active_front = _neighbor_front(i);
                Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : domain_idx;
                Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : domain_idx;
                Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : domain_idx;
                Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : domain_idx;
                
                // Interpolate u to y-face (uy)
                _velocity_x_at_y(domain_idx) = 0.25 * (_velocity_x(domain_idx) + _velocity_x(i_left) +
                                                       _velocity_x(i_front) + 
                                                       (active_front >= 0 && _neighbor_left(active_front) >= 0 ? 
                                                        _velocity_x(_active_to_domain(_neighbor_left(active_front))) : 
                                                        _velocity_x(i_front)));
                
                // Interpolate v to x-face (vx)
                _velocity_y_at_x(domain_idx) = 0.25 * (_velocity_y(domain_idx) + _velocity_y(i_back) +
                                                       _velocity_y(i_right) +
                                                       (active_right >= 0 && _neighbor_back(active_right) >= 0 ?
                                                        _velocity_y(_active_to_domain(_neighbor_back(active_right))) :
                                                        _velocity_y(i_right)));
            });
    }
    
    // ========================================================================
    // APPLY CFL LIMITER
    // ========================================================================
    void apply_cfl_limiter() {
        auto _pressure = state.pressure;
        auto _bottom = state.bottom;
        auto _depth = state.depth;
        auto _cfl_active = state.cfl_active;
        
        auto _neighbor_right = active_mesh.neighbor_right;
        auto _neighbor_left = active_mesh.neighbor_left;
        auto _neighbor_front = active_mesh.neighbor_front;
        auto _neighbor_back = active_mesh.neighbor_back;
        auto _active_to_domain = active_mesh.active_to_domain;
        
        Kokkos::parallel_for(RangePolicy(0, active_mesh.num_active),
            [=] KOKKOS_INLINE_FUNCTION (const Ordinal i) {
                Ordinal domain_idx = _active_to_domain(i);
                
                Scalar diff = _pressure(domain_idx) - _bottom(domain_idx);
                
                // Check for isolated wet cells
                if (_depth(domain_idx) <= 0.0 && diff > 0.0) {
                    bool has_wet_neighbor = false;
                    
                    Ordinal active_right = _neighbor_right(i);
                    Ordinal active_left = _neighbor_left(i);
                    Ordinal active_front = _neighbor_front(i);
                    Ordinal active_back = _neighbor_back(i);
                    Ordinal i_right = (active_right >= 0) ? _active_to_domain(active_right) : -1;
                    Ordinal i_left = (active_left >= 0) ? _active_to_domain(active_left) : -1;
                    Ordinal i_front = (active_front >= 0) ? _active_to_domain(active_front) : -1;
                    Ordinal i_back = (active_back >= 0) ? _active_to_domain(active_back) : -1;
                    
                    if (i_right >= 0 && _depth(i_right) > 0.0) has_wet_neighbor = true;
                    if (i_left >= 0 && _depth(i_left) > 0.0) has_wet_neighbor = true;
                    if (i_front >= 0 && _depth(i_front) > 0.0) has_wet_neighbor = true;
                    if (i_back >= 0 && _depth(i_back) > 0.0) has_wet_neighbor = true;
                    
                    if (!has_wet_neighbor) {
                        _pressure(domain_idx) = _bottom(domain_idx);
                        _cfl_active(domain_idx) = 1;
                    }
                }
                
                // Remove small depth
                if (diff > 0.0 && diff < min_depth) {
                    _pressure(domain_idx) = _bottom(domain_idx);
                }
            });
    }
};

#endif // FREHG_SHALLOW_WATER_SOLVER_HPP

