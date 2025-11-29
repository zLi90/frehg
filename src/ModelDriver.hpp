#ifndef FREHG_MODEL_DRIVER_HPP
#define FREHG_MODEL_DRIVER_HPP

#include "define.hpp"
#include <mpi.h>
#include <Kokkos_Core.hpp> // Required for Kokkos::Timer
#include <string>
#include <iostream>
#include <memory> 

// Forward declarations of future modules
// class Mesh;
// class SurfaceFlow;
// class SubsurfaceFlow;
// class Transport;

namespace Frehg {

class ModelDriver {
private:
    // MPI Communicator for this simulation instance
    MPI_Comm comm_;
    int rank_;
    
    // Performance timer
    Kokkos::Timer timer_;
    
    // Future pointers to modules
    // std::unique_ptr<Mesh> mesh_;
    // std::unique_ptr<SubsurfaceFlow> gw_model_;
    // std::unique_ptr<SurfaceFlow> sw_model_;

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
        if (rank_ == 0) {
            std::cout << "------------------------------------------------------------\n";
            std::cout << "[Driver] Phase 1: Initialization\n";
            std::cout << "[Driver] Reading configuration from: " << input_dir << "\n";
            std::cout << "[Driver] Writing results to        : " << output_dir << "\n";
            std::cout << "------------------------------------------------------------\n";
        }

        // 1. TODO: Parse Input File (e.g., input_dir/frehg_config.txt)
        //    Load time step size, total simulation time, physics flags.

        // 2. TODO: Initialize Mesh
        //    mesh_ = std::make_unique<Mesh>(input_dir + "/mesh_file.exo");
        
        // 3. TODO: Initialize Physics Modules
        //    gw_model_ = std::make_unique<SubsurfaceFlow>(mesh_.get(), parameters);
        
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
        
        // Placeholder time stepping logic
        Scalar current_time = 0.0;
        Scalar end_time = 100.0; // This would come from input file
        Scalar dt = 1.0;         // This would come from input file
        int step = 0;

        while (current_time < end_time) {
            // 1. Solve Surface (Shallow Water)
            // sw_model_->solve(dt);
            
            // 2. Exchange Fluxes (Coupling)
            // exchange_fluxes();
            
            // 3. Solve Subsurface (Richards)
            // gw_model_->solve(dt);
            
            // 4. Solve Transport (Salinity)
            // transport_->solve(dt);
            
            // 5. I/O (Write output files at specified intervals)
            // if (step % output_interval == 0) write_output(step);
            
            // Update time
            current_time += dt;
            step++;
            
            // Simple progress log
            if (rank_ == 0 && step % 10 == 0) {
                std::cout << "Step: " << step << " | Time: " << current_time << "s\r" << std::flush;
            }
        }

        // Barrier to ensure all ranks finish before printing time
        MPI_Barrier(comm_);
        
        double total_time = timer_.seconds();
        if (rank_ == 0) {
            std::cout << "\n[Driver] Finished. Total wall time: " << total_time << " s" << std::endl;
        }
    }
};

} // namespace Frehg

#endif // FREHG_MODEL_DRIVER_HPP