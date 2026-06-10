#include "core/Simulation.hpp"
#include "driver/SimulationDriver.hpp"

#include <exception>
#include <iostream>
#include <string>

#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_PETSC
#include <petsc.h>
#endif

namespace {

bool isHelpRequest(int argc, char** argv)
{
    if (argc != 2 || argv[1] == nullptr) {
        return false;
    }

    const std::string argument = argv[1];
    return argument == "--help" || argument == "-h";
}

}  // namespace

int main(int argc, char** argv)
{
    if (isHelpRequest(argc, argv)) {
        std::cout << frehg2::helpText(argv[0]);
        return 0;
    }

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, nullptr, nullptr);
#elif defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif

#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#endif

    int exit_code = 0;
    try {
        const auto config = frehg2::makeConfigFromCommandLine(argc, argv);
        exit_code = frehg2::SimulationDriver(config.input_path).run();
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n\n";
        std::cerr << frehg2::helpText(argv[0]);
        exit_code = 1;
    }

#ifdef USE_KOKKOS
    Kokkos::finalize();
#endif

#ifdef USE_PETSC
    PetscFinalize();
#elif defined(USE_MPI)
    MPI_Finalize();
#endif

    return exit_code;
}
