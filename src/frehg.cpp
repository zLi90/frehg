/*
 * File: frehg.cpp
 * Description: Main entry point for the Frehg coupled model.
 * Usage: mpirun -n <N> ./frehg <input_dir> <output_dir> <num_threads>
 */

#include <iostream>
#include <string>
#include <cstdlib> // for std::stoi
#include <mpi.h>
#include <Kokkos_Core.hpp>

#include "define.hpp"
#include "ModelDriver.hpp"

int main(int argc, char* argv[]) {
    // 1. Initialize MPI
    // ------------------------------------------------------------------------
    int mpi_provided;
    // MPI_THREAD_FUNNELED implies only the main thread will make MPI calls, 
    // which is standard for Hybrid MPI+Kokkos.
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided);

    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // 2. Parse Runtime Arguments
    // ------------------------------------------------------------------------
    if (argc != 4) {
        if (mpi_rank == 0) {
            std::cerr << "Error: Invalid number of arguments.\n";
            std::cerr << "Usage: mpirun -n <NP> ./frehg <input_dir> <output_dir> <num_threads>\n";
        }
        MPI_Finalize();
        return 1;
    }

    std::string input_dir = argv[1];
    std::string output_dir = argv[2];
    int num_threads = std::stoi(argv[3]);

    // 3. Initialize Kokkos with Explicit Thread Count
    // ------------------------------------------------------------------------
    // Kokkos::InitArguments allows us to pass runtime configuration 
    // without relying on environment variables or --kokkos-flags.
    Kokkos::InitArguments args;
    
    // Set the number of threads for the Host execution space (e.g., OpenMP).
    // Note: If running on GPU, this controls the CPU-side threads, 
    // while the GPU uses its own massive parallelism.
    args.num_threads = num_threads;
    
    // We pass 'true' for warnings to catch any configuration issues (e.g., oversubscription)
    Kokkos::initialize(args);

    // 4. Simulation Scope (RAII)
    // ------------------------------------------------------------------------
    // The brace block ensures 'driver' is destroyed before Kokkos::finalize()
    {
        if (mpi_rank == 0) {
            std::cout << "============================================================\n";
            std::cout << "   Frehg: Fine-Resolution Environmental Hydrodynamic & GW   \n";
            std::cout << "============================================================\n";
            std::cout << "MPI Size         : " << mpi_size << "\n";
            std::cout << "Requested Threads: " << num_threads << "\n";
            std::cout << "Kokkos ExecSpace : " << typeid(Kokkos::DefaultExecutionSpace).name() << "\n";
            std::cout << "Input Directory  : " << input_dir << "\n";
            std::cout << "Output Directory : " << output_dir << "\n";
            std::cout << "============================================================\n" << std::endl;
        }

        try {
            // Instantiate the Model Driver
            Frehg::ModelDriver driver(MPI_COMM_WORLD);

            // Phase 1: Setup
            // We pass the parsed directories to the driver
            driver.initialize(input_dir, output_dir);

            // Phase 2: Execution
            driver.run();

            if (mpi_rank == 0) {
                std::cout << "\nSimulation completed successfully." << std::endl;
            }

        } catch (const std::exception& e) {
            std::cerr << "[Rank " << mpi_rank << "] Exception: " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    // 5. Finalize Resources
    // ------------------------------------------------------------------------
    Kokkos::finalize();
    MPI_Finalize();

    return 0;
}