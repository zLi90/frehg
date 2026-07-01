// Frehg2 driver entry point (P1.4 / P7).
//
// Brings up and tears down the three frameworks in the required order
// (Kokkos -> MPI -> PETSc on init; PETSc -> Kokkos -> MPI on shutdown), parses the command
// line, and runs the simulation through the Orchestrator — the ONLY production path from
// main() to the physics. --help / --version are handled before any framework init so they
// work without MPI.
#include <petscsys.h>

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <cstring>
#include <exception>
#include <iostream>
#include <string>

#include "driver/Orchestrator.hpp"
#include "io/Config.hpp"

namespace {

constexpr const char* kProgramName = "frehg2";
constexpr const char* kVersion = "0.1.0";

void printUsage(std::ostream& os) {
  os << "Usage: " << kProgramName << " [options] [config.yaml]\n"
     << "\n"
     << "  Surface water / groundwater / solute transport coupled numerical model.\n"
     << "\n"
     << "Options:\n"
     << "  -h, --help                 Show this help message and exit\n"
     << "  -v, --version              Show version information and exit\n"
     << "  --restart FILE             Restart from an HDF5 checkpoint file\n"
     << "  --restart-time T           Simulation time (s) stored in the checkpoint\n"
     << "\n"
     << "PETSc options (e.g. -ksp_type, -mat_type) may be passed through and are\n"
     << "consumed by PetscInitialize.\n";
}

// Value of a "--flag value" option (empty string if absent).
std::string optionValue(int argc, char** argv, const char* flag) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (std::strcmp(argv[i], flag) == 0) return std::string(argv[i + 1]);
  }
  return std::string();
}

bool hasFlag(int argc, char** argv, const char* shortFlag, const char* longFlag) {
  for (int i = 1; i < argc; ++i) {
    if ((shortFlag && std::strcmp(argv[i], shortFlag) == 0) ||
        (longFlag && std::strcmp(argv[i], longFlag) == 0)) {
      return true;
    }
  }
  return false;
}

// First non-flag argument (that is not the value of a "--flag value" option) is the YAML
// config path.
std::string findConfigPath(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    if (i > 0 && (std::strcmp(argv[i - 1], "--restart") == 0 ||
                  std::strcmp(argv[i - 1], "--restart-time") == 0)) {
      continue;  // this token is an option value, not the config path
    }
    if (argv[i][0] != '-') {
      return std::string(argv[i]);
    }
  }
  return std::string();
}

bool fileExists(const std::string& path) {
  if (path.empty()) {
    return false;
  }
  FILE* f = std::fopen(path.c_str(), "r");
  if (f != nullptr) {
    std::fclose(f);
    return true;
  }
  return false;
}

}  // namespace

int main(int argc, char** argv) {
  // --help / --version before any framework initialization.
  if (hasFlag(argc, argv, "-h", "--help")) {
    printUsage(std::cout);
    return 0;
  }
  if (hasFlag(argc, argv, "-v", "--version")) {
    std::cout << kProgramName << " version " << kVersion << '\n';
    return 0;
  }

  // Framework init order: Kokkos first (so device views exist before PETSc queries),
  // then MPI, then PETSc.
  Kokkos::initialize(argc, argv);

  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cerr << kProgramName << ": MPI_Init failed\n";
    Kokkos::finalize();
    return 1;
  }

  PetscErrorCode petsc_ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
  if (petsc_ierr != 0) {
    std::cerr << kProgramName << ": PetscInitialize failed (code " << petsc_ierr << ")\n";
    MPI_Finalize();
    Kokkos::finalize();
    return 1;
  }

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int exit_code = 0;
  const std::string config_path = findConfigPath(argc, argv);

  if (config_path.empty()) {
    if (rank == 0) {
      std::cerr << kProgramName << ": error: no config file given\n";
      printUsage(std::cerr);
    }
    exit_code = 1;
  } else if (!fileExists(config_path)) {
    if (rank == 0) {
      std::cerr << kProgramName << ": error: config file not found: " << config_path << '\n';
    }
    exit_code = 1;
  } else {
    try {
      frehg2::Config cfg = frehg2::Config::fromFile(config_path);
      frehg2::Orchestrator orch;
      orch.initialize(cfg);
      const std::string restart_file = optionValue(argc, argv, "--restart");
      if (!restart_file.empty()) {
        const std::string rt = optionValue(argc, argv, "--restart-time");
        const double restart_time = rt.empty() ? 0.0 : std::stod(rt);
        orch.restart(restart_file, restart_time);
      } else {
        orch.run();
      }
      if (rank == 0) {
        std::cout << kProgramName << ": run '" << config_path << "' complete ("
                  << orch.stepCount() << " steps; summary -> " << orch.summaryPath() << ")\n";
      }
    } catch (const std::exception& e) {
      if (rank == 0) std::cerr << kProgramName << ": error: " << e.what() << '\n';
      exit_code = 1;
    }
  }

  // Shutdown: PETSc, then Kokkos, then MPI (PETSc must finalize before MPI).
  PetscFinalize();
  Kokkos::finalize();
  MPI_Finalize();
  return exit_code;
}
