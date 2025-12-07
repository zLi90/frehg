/*
 * File: frehg.cpp
 * Description: Main entry point for the Frehg coupled model.
 * Usage: ./frehg <input_dir> <output_dir> [num_threads]
 */

#include <iostream>
#include <string>
#include <cstdlib> // for std::stoi
#include <Kokkos_Core.hpp>

#include "define.hpp"
#include "ModelDriver.hpp"

int main(int argc, char* argv[]) {
    // 1. Parse Runtime Arguments
    // ------------------------------------------------------------------------
    if (argc < 3 || argc > 4) {
        std::cerr << "Error: Invalid number of arguments.\n";
        std::cerr << "Usage: ./frehg <input_dir> <output_dir> [num_threads]\n";
        std::cerr << "  num_threads: Optional, defaults to 1\n";
        return 1;
    }

    std::string input_dir = argv[1];
    std::string output_dir = argv[2];
    int num_threads = (argc >= 4) ? std::stoi(argv[3]) : 1;

    // 2. Initialize Kokkos with Explicit Thread Count
    // ------------------------------------------------------------------------
    // Set the number of threads for the Host execution space (e.g., OpenMP).
    Kokkos::InitializationSettings settings;
    settings.set_num_threads(num_threads);
    
    Kokkos::initialize(settings);

    // 3. Simulation Scope (RAII)
    // ------------------------------------------------------------------------
    // The brace block ensures 'driver' is destroyed before Kokkos::finalize()
    {
        std::cout << "============================================================\n";
        std::cout << "   Frehg: Fine-Resolution Environmental Hydrodynamic & GW   \n";
        std::cout << "============================================================\n";
        std::cout << "Kokkos Threads   : " << num_threads << "\n";
        std::cout << "Kokkos ExecSpace : " << typeid(Kokkos::DefaultExecutionSpace).name() << "\n";
        std::cout << "Input Directory  : " << input_dir << "\n";
        std::cout << "Output Directory : " << output_dir << "\n";
        std::cout << "============================================================\n" << std::endl;

        try {
            // Instantiate the Model Driver
            Frehg::ModelDriver driver;

            // Phase 1: Setup
            // We pass the parsed directories to the driver
            driver.initialize(input_dir, output_dir);

            // Phase 2: Execution
            driver.run();

            std::cout << "\nSimulation completed successfully." << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            Kokkos::finalize();
            return 1;
        }
    }

    // 4. Finalize Resources
    // ------------------------------------------------------------------------
    Kokkos::finalize();

    return 0;
}
