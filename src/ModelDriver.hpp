#ifndef FREHG_MODEL_DRIVER_HPP
#define FREHG_MODEL_DRIVER_HPP

#include "define.hpp"
#include "Domain.hpp"
#include "StateVariables.hpp"
#include "ActiveCellMesh.hpp"
#include "ShallowWaterSolver.hpp"
#include "GroundwaterSolver.hpp"
#include "ScalarTransportSolver.hpp"
#include "Initializer.hpp"
#include <mpi.h>
#include <Kokkos_Core.hpp> // Required for Kokkos::Timer
#include <string>
#include <iostream>
#include <memory> 
#include <vector>
#include <cmath>

namespace Frehg {

class ModelDriver {
private:
    // MPI Communicator for this simulation instance
    MPI_Comm comm_;
    int rank_;
    
    // Performance timer
    Kokkos::Timer timer_;
    
    // Note: Domains, state variables, meshes, and solvers are owned by initializer_
    // Access them via initializer_->get_*() methods
    
    // --- Simulation Parameters ---
    bool sim_shallowwater_;
    bool sim_groundwater_;
    bool sync_coupling_;          // Synchronous coupling (same dt) vs async
    int n_scalar_;                // Number of scalar species
    Scalar dt_;                   // Surface water time step
    Scalar dtg_;                  // Groundwater time step (may differ)
    Scalar dt_min_;               // Minimum groundwater time step for async
    Scalar end_time_;             // Total simulation time
    Scalar dt_out_;               // Output interval
    int output_interval_;          // Output every N steps
    
    // --- Output Directory ---
    std::string output_dir_;
    
    // --- Initializer ---
    std::unique_ptr<Initializer> initializer_;

public:
    // ========================================================================
    // Constructor
    // ========================================================================
    ModelDriver(MPI_Comm comm) : comm_(comm) {
        MPI_Comm_rank(comm_, &rank_);
    }

    // ========================================================================
    // Destructor
    // ========================================================================
    ~ModelDriver() = default;

    // ========================================================================
    // Initialization Phase
    // ========================================================================
    void initialize(const std::string& input_dir, const std::string& output_dir) {
        output_dir_ = output_dir;
        
        if (rank_ == 0) {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "[Driver] Phase 1: Initialization\n";
            std::cout << "[Driver] Reading configuration from: " << input_dir << "\n";
            std::cout << "[Driver] Writing results to        : " << output_dir << "\n";
            std::cout << "------------------------------------------------------------\n";
        }

        // Create initializer and perform all initialization steps
        initializer_ = std::make_unique<Initializer>(input_dir, output_dir);
        initializer_->initialize_all();
        
        // Extract configuration from initializer
        const auto& config = initializer_->get_config();
        sim_shallowwater_ = config.sim_shallowwater;
        sim_groundwater_ = config.sim_groundwater;
        sync_coupling_ = config.sync_coupling;
        n_scalar_ = config.n_scalar;
        dt_ = config.dt;
        dtg_ = config.dt;
        dt_min_ = config.dt * 0.1;  // Default: 10% of main dt
        end_time_ = config.Tend;
        dt_out_ = config.dt_out;
        output_interval_ = static_cast<int>(config.dt_out / config.dt);
        
        // Objects are owned by initializer_, accessed via getters
        
        if (rank_ == 0) {
            std::cout << "[Driver] Initialization complete.\n" << std::endl;
        }
    }

    // ========================================================================
    // Execution Phase (Time Loop)
    // ========================================================================
    void run() {
        if (rank_ == 0) {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "[Driver] Phase 2: Simulation Loop\n";
            std::cout << "------------------------------------------------------------\n";
        }
        
        // Reset timer to measure pure computation time (excluding setup)
        timer_.reset();
        
        // Initialize time stepping variables
        Scalar current_time = 0.0;
        Scalar t_subsurface = 0.0;  // Subsurface time (for async coupling)
        Scalar last_save = 0.0;
        int step = 0;
        
        // Set groundwater time step
        if (!sim_groundwater_ || sync_coupling_) {
            dtg_ = dt_;
        } else {
            dtg_ = dt_min_;
        }
        
        // Write initial conditions
        if (rank_ == 0) {
            std::cout << "[Driver] Writing initial conditions...\n";
        }
        // TODO: write_output(0, current_time);
        
        if (rank_ == 0) {
            std::cout << "[Driver] Beginning time loop...\n";
        }

        while (current_time < end_time_) {
            // Update time
            current_time += dt_;
            if (!sim_groundwater_ || sync_coupling_) {
                t_subsurface = current_time;
            }
            step++;
            
            if (rank_ == 0 && step % 10 == 0) {
                std::cout << "Step: " << step << " | Time: " << current_time 
                         << "s / " << end_time_ << "s\r" << std::flush;
            }
            
            // 1. Get boundary conditions at current time
            // TODO: get_current_bc(current_time);
            
            // 2. Get evaporation/rainfall at current time
            // TODO: get_evaprain(current_time);
            
            // 3. Solve Shallow Water (if enabled)
            if (sim_shallowwater_ && initializer_) {
                auto* solver = initializer_->get_sw_solver();
                if (solver) solver->solve(current_time);
            }
            
            // 4. Solve Groundwater (if enabled)
            if (sim_groundwater_ && initializer_) {
                // Compute seepage rate for top boundary (if coupling)
                // TODO: compute_seepage_rate();
                
                auto* solver = initializer_->get_gw_solver();
                if (solver) {
                    if (sync_coupling_) {
                        // Synchronous coupling: same time step
                        solver->solve(current_time);
                    } else {
                        // Asynchronous coupling: multiple groundwater steps per surface step
                        while (t_subsurface + dtg_ <= current_time) {
                            t_subsurface += dtg_;
                            solver->solve(t_subsurface);
                            if (rank_ == 0) {
                                std::cout << "   >>> Subsurface executed with dt = " << dtg_ << "\n";
                            }
                        }
                    }
                }
            }
            
            // 5. Update Shallow Water Velocity (after both solvers)
            if (sim_shallowwater_ && initializer_) {
                auto* solver = initializer_->get_sw_solver();
                if (solver) solver->update_velocity();
            }
            
            // 6. Check CFL number (for surface water)
            // TODO: check_cfl();
            
            // 7. Solve Scalar Transport (if enabled)
            if (n_scalar_ > 0 && initializer_) {
                // Surface water scalar transport
                if (sim_shallowwater_) {
                    const auto& solvers = initializer_->get_sw_scalar_solvers();
                    for (size_t kk = 0; kk < solvers.size() && kk < static_cast<size_t>(n_scalar_); ++kk) {
                        if (solvers[kk]) {
                            solvers[kk]->solve(current_time);
                        }
                    }
                }
                
                // Groundwater scalar transport
                if (sim_groundwater_) {
                    const auto& solvers = initializer_->get_gw_scalar_solvers();
                    for (size_t kk = 0; kk < solvers.size() && kk < static_cast<size_t>(n_scalar_); ++kk) {
                        if (solvers[kk]) {
                            solvers[kk]->solve(current_time);
                        }
                    }
                }
            }
            
            // 8. Reset rainfall (if needed)
            // TODO: reset_rainfall();
            
            // 9. Write output at specified intervals
            if (std::abs(current_time - last_save - dt_out_) <= dt_) {
                int t_save = static_cast<int>(std::round(current_time / dt_out_)) * static_cast<int>(dt_out_);
                last_save = t_save;
                // TODO: write_output(step, current_time);
            }
            
            // 10. Write monitored variables
            // TODO: write_monitor_output(step, current_time);
            
            // 11. Report time step completion
            if (rank_ == 0 && step % 100 == 0) {
                std::cout << "\n[Driver] Step " << step << " (" << current_time 
                         << "s of " << end_time_ << "s) completed, dt = " << dt_ << "\n";
            }
        }

        // Barrier to ensure all ranks finish before printing time
        MPI_Barrier(comm_);
        
        double total_time = timer_.seconds();
        if (rank_ == 0) {
            std::cout << "\n[Driver] Simulation completed!\n";
            std::cout << "[Driver] Total steps: " << step << "\n";
            std::cout << "[Driver] Total wall time: " << total_time << " s\n";
            std::cout << "[Driver] Average time per step: " << total_time / step << " s\n";
        }
    }
    
private:
    // ========================================================================
    // Helper Functions (Placeholders for future implementation)
    // ========================================================================
    
    // Get boundary conditions at current time
    void get_current_bc(Scalar t_current) {
        // TODO: Implement boundary condition interpolation
        // - Wind conditions
        // - Tidal boundary conditions
        // - Inflow boundary conditions
        // - Scalar boundary conditions
    }
    
    // Get evaporation/rainfall at current time
    void get_evaprain(Scalar t_current) {
        // TODO: Implement evaporation/rainfall update
        // - Interpolate from time series
        // - Compute aerodynamic evaporation if needed
        // - Update qtop for groundwater
    }
    
    // Compute seepage rate for surface-subsurface coupling
    void compute_seepage_rate() {
        // TODO: Compute seepage flux at top boundary
        // if (sim_shallowwater_ && sim_groundwater_) {
        //     // Compute rss based on surface elevation change
        // }
    }
    
    // Check CFL number
    void check_cfl() {
        // TODO: Compute and check CFL numbers
        // - Get max CFL from surface water solver
        // - Optionally adjust time step if CFL > 1.0
    }
    
    // Write output files
    void write_output(int step, Scalar time) {
        // TODO: Write output files
        // - Surface water fields (eta, depth, velocity, scalar)
        // - Groundwater fields (head, water content, scalar)
        // - Use output_dir_ for file paths
    }
    
    // Write monitored variables
    void write_monitor_output(int step, Scalar time) {
        // TODO: Write time series for monitored locations
    }
    
    // Reset rainfall
    void reset_rainfall() {
        // TODO: Reset rainfall accumulation if needed
    }
};

} // namespace Frehg

#endif // FREHG_MODEL_DRIVER_HPP