cmake_minimum_required(VERSION 3.20)

project(Frehg2
  VERSION 0.1.0
  DESCRIPTION "Surface water / groundwater / solute transport coupled numerical model"
  LANGUAGES C CXX)

# ---------------------------------------------------------------------------------------
# Language standard (C++20: baseline for modern Kokkos; matches the local Kokkos build
# which was configured with Kokkos_CXX_STANDARD=20 and g++-15).
# ---------------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------------------
option(BUILD_TESTING "Build the test suite and register CTest tests" ON)
option(FREHG2_STRICT_WARNINGS "Treat compiler warnings as errors (-Wall -Wextra -Werror)" ON)
option(FREHG2_ENABLE_CUDA "Enable the CUDA GPU backend (Linux/NVIDIA only)" OFF)
option(FREHG2_USE_GPU_MAT "Default to MATAIJKOKKOS/VECKOKKOS PETSc types (override with -mat_type)" OFF)
# P20: the legacy frehg C/MPI/LASPack code is DEPRECATED and NOT compiled by default. It is kept
# under legacy/frehg/ for archival/reproducibility only (see legacy/frehg/DEPRECATED.md). Set ON to
# build it via its own makefile (archival opt-in); Frehg2 itself never links or includes it.
option(FREHG2_USE_LEGACY "Build the archived legacy frehg code via its own makefile (opt-in)" OFF)

# Dependency prefix for the local toolchain (PETSc, Kokkos, HDF5, yaml-cpp, MPICH).
set(FREHG2_LOCAL_PREFIX "/Users/zhili/Codes/local" CACHE PATH "Local dependency install prefix")

# macOS cannot run CUDA: fail fast and loud if the GPU backend is requested here.
if(APPLE AND FREHG2_ENABLE_CUDA)
  message(FATAL_ERROR
    "FREHG2_ENABLE_CUDA=ON is not supported on macOS (no CUDA). "
    "Build the production CPU/OpenMP model here and defer GPU execution to a "
    "Linux/NVIDIA machine (see docs/gpu_validation_policy.md).")
endif()

# Guard against an accidental compiler switch: the plan/.cursorrules pin g++-15 (which is
# also what Kokkos was built with). Warn loudly if a different compiler is in use.
if(NOT CMAKE_CXX_COMPILER MATCHES "g\\+\\+-15$")
  message(WARNING
    "CMAKE_CXX_COMPILER is '${CMAKE_CXX_COMPILER}', expected g++-15 (the compiler the "
    "local Kokkos was built with). Configure with -DCMAKE_CXX_COMPILER=g++-15.")
endif()

# ---------------------------------------------------------------------------------------
# Dependency discovery
# ---------------------------------------------------------------------------------------
# MPI: pin to the local MPICH that PETSc was built against (NOT Homebrew Open MPI). See
# docs/plan_reconciliation.md (dependency layout, MPICH risk).
if(NOT MPI_CXX_COMPILER AND EXISTS "${FREHG2_LOCAL_PREFIX}/bin/mpicxx")
  set(MPI_CXX_COMPILER "${FREHG2_LOCAL_PREFIX}/bin/mpicxx" CACHE FILEPATH "MPI C++ wrapper (MPICH)")
endif()
if(NOT MPI_C_COMPILER AND EXISTS "${FREHG2_LOCAL_PREFIX}/bin/mpicc")
  set(MPI_C_COMPILER "${FREHG2_LOCAL_PREFIX}/bin/mpicc" CACHE FILEPATH "MPI C wrapper (MPICH)")
endif()
if(NOT MPIEXEC_EXECUTABLE AND EXISTS "${FREHG2_LOCAL_PREFIX}/bin/mpiexec")
  set(MPIEXEC_EXECUTABLE "${FREHG2_LOCAL_PREFIX}/bin/mpiexec" CACHE FILEPATH "MPICH mpiexec")
endif()
find_package(MPI REQUIRED)

# PETSc: via cmake/FindPETSc.cmake (pkg-config wrapper; no upstream CMake config).
if(NOT PETSc_DIR)
  set(PETSc_DIR "${FREHG2_LOCAL_PREFIX}" CACHE PATH "PETSc install prefix")
endif()
find_package(PETSc REQUIRED)

# Kokkos: upstream CMake config package.
if(NOT Kokkos_ROOT)
  set(Kokkos_ROOT "${FREHG2_LOCAL_PREFIX}" CACHE PATH "Kokkos install prefix")
endif()
find_package(Kokkos REQUIRED)

# HDF5: C API only (no H5Cpp.h locally; hdf5.pc is stale). Discover via HDF5_ROOT.
if(NOT HDF5_ROOT)
  set(HDF5_ROOT "${FREHG2_LOCAL_PREFIX}" CACHE PATH "HDF5 install prefix")
endif()
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 REQUIRED COMPONENTS C)

# yaml-cpp: upstream CMake config package.
if(NOT yaml-cpp_ROOT)
  set(yaml-cpp_ROOT "${FREHG2_LOCAL_PREFIX}" CACHE PATH "yaml-cpp install prefix")
endif()
find_package(yaml-cpp REQUIRED)

# ---------------------------------------------------------------------------------------
# Discovery report (configure-time verification)
# ---------------------------------------------------------------------------------------
message(STATUS "==================== Frehg2 ${PROJECT_VERSION} configuration ====================")
message(STATUS "  C++ standard      : C++${CMAKE_CXX_STANDARD}")
message(STATUS "  CXX compiler      : ${CMAKE_CXX_COMPILER}")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  MPI (CXX) wrapper : ${MPI_CXX_COMPILER}")
message(STATUS "  MPI include       : ${MPI_CXX_INCLUDE_DIRS}")
message(STATUS "  PETSc ${PETSC_VERSION} : ${PETSC_LIBRARY}")
message(STATUS "  PETSc include     : ${PETSC_INCLUDE_DIR}")
message(STATUS "  Kokkos ${Kokkos_VERSION} : ${Kokkos_DIR}")
message(STATUS "  Kokkos devices    : ${Kokkos_DEVICES}")
message(STATUS "  HDF5 ${HDF5_VERSION} : ${HDF5_C_LIBRARIES}")
message(STATUS "  yaml-cpp          : ${yaml-cpp_DIR}")
message(STATUS "  STRICT_WARNINGS   : ${FREHG2_STRICT_WARNINGS}")
message(STATUS "  ENABLE_CUDA       : ${FREHG2_ENABLE_CUDA} (CUDA disabled on macOS)")
message(STATUS "  USE_GPU_MAT       : ${FREHG2_USE_GPU_MAT}")
message(STATUS "  LEGACY (frehg)    : ${FREHG2_USE_LEGACY} (DEPRECATED; not compiled unless ON)")
message(STATUS "================================================================================")

# ---------------------------------------------------------------------------------------
# Legacy archival build (P20 — Legacy Deprecation & Removal)
# ---------------------------------------------------------------------------------------
# The legacy serial C / MPI / LASPack `frehg` code served as the numerical ground truth through
# P0-P19. Now that the unified validation harness (P19) is the canonical test, the legacy code is
# DEPRECATED: it is NOT part of the Frehg2 build, is never linked or #included by any Frehg2 target,
# and is compiled ONLY when a user explicitly opts in for archival/reproducibility. It builds with
# its OWN historical build system (legacy/frehg/makefile, mpicc + LASPack), not CMake.
if(FREHG2_USE_LEGACY)
  message(WARNING
    "FREHG2_USE_LEGACY=ON: building the DEPRECATED legacy frehg code via legacy/frehg/makefile. "
    "This is an archival / reproducibility path only -- Frehg2 is the supported model. "
    "See legacy/frehg/DEPRECATED.md.")
  find_program(MAKE_EXECUTABLE NAMES gmake make REQUIRED)
  add_custom_target(legacy_frehg ALL
    COMMAND ${MAKE_EXECUTABLE} -C "${CMAKE_SOURCE_DIR}/legacy/frehg"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/legacy/frehg"
    COMMENT "Building archived legacy frehg (DEPRECATED) via its own makefile"
    VERBATIM)
endif()

# ---------------------------------------------------------------------------------------
# Code-quality helpers (CQ1, CQ6)
# ---------------------------------------------------------------------------------------
# CQ1: strict warnings applied uniformly to every Frehg2 target.
function(apply_strict_warnings target)
  if(FREHG2_STRICT_WARNINGS)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Werror -Wno-unused-parameter)
  endif()
endfunction()

# CQ6: forbid sequential-only PETSc patterns leaking into production code.
add_custom_target(check_no_seqaij
  COMMAND ${CMAKE_COMMAND} -E echo "Checking for sequential PETSc patterns in src/ ..."
  COMMAND bash -c "! grep -rn 'MatCreateSeqAIJ\\|MatCreateSeqBAIJ' ${CMAKE_SOURCE_DIR}/src/"
  COMMAND bash -c "! grep -rn 'PETSC_COMM_SELF' ${CMAKE_SOURCE_DIR}/src/"
  COMMENT "Verify no sequential PETSc patterns (MatCreateSeqAIJ / PETSC_COMM_SELF) in src/"
  VERBATIM)

# P2.6.3: enforce the solver seam. (1) MatSetValues only inside the PETSc backend bridge;
# (2) no solver-library types leak into physics directories.
add_custom_target(check_solver_seam
  COMMAND ${CMAKE_COMMAND} -E echo "Checking the solver seam (no solver types in physics)..."
  COMMAND bash -c "! grep -rn 'MatSetValues' ${CMAKE_SOURCE_DIR}/src | grep -v 'linear/backends/'"
  COMMAND bash -c "! grep -rEn 'Mat |Vec |KSP |DM |DMDA|Tpetra' ${CMAKE_SOURCE_DIR}/src/swe ${CMAKE_SOURCE_DIR}/src/re ${CMAKE_SOURCE_DIR}/src/solute ${CMAKE_SOURCE_DIR}/src/coupling 2>/dev/null"
  COMMENT "Verify physics code uses SparseSystem/LinearSolver, not solver-library types"
  VERBATIM)

# P10 acceptance: the GPU-native assembly TU must use MatSetValuesCOO ONLY — no host-staging
# MatSetValues call may appear in the device code path. (MatSetValuesCOO is allowed; the grep
# excludes it explicitly.) This is the static check behind "MatSetValues absent from the GPU
# code path".
add_custom_target(check_gpu_coo_assembly
  COMMAND ${CMAKE_COMMAND} -E echo "Checking the GPU assembly path uses MatSetValuesCOO only..."
  COMMAND bash -c "! grep -nE 'MatSetValues\\(' ${CMAKE_SOURCE_DIR}/src/linear/backends/KokkosPetscBridgeDevice.cpp"
  COMMAND bash -c "grep -q 'MatSetValuesCOO' ${CMAKE_SOURCE_DIR}/src/linear/backends/KokkosPetscBridgeDevice.cpp"
  COMMENT "Verify the device assembly path is MatSetValuesCOO-only (no host MatSetValues)"
  VERBATIM)

# Convenience interface target for the optional GPU-matrix default.
if(FREHG2_USE_GPU_MAT)
  add_compile_definitions(FREHG2_USE_GPU_MAT=1)
endif()

# P10: surface the CUDA build intent as a compile definition (the actual assembly-path switch
# keys on the Kokkos device macros via GpuAssembly.hpp, so this is informational/diagnostic).
if(FREHG2_ENABLE_CUDA)
  add_compile_definitions(FREHG2_ENABLE_CUDA=1)
endif()

# ---------------------------------------------------------------------------------------
# Build provenance (P3.2.1): bake git SHA / build info into a generated header.
# ---------------------------------------------------------------------------------------
find_package(Git QUIET)
set(FREHG2_GIT_SHA "unknown")
set(FREHG2_GIT_DIRTY "false")
if(Git_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short=12 HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE _git_sha OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET RESULT_VARIABLE _git_rc)
  if(_git_rc EQUAL 0 AND _git_sha)
    set(FREHG2_GIT_SHA "${_git_sha}")
    execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE _git_status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(_git_status)
      set(FREHG2_GIT_DIRTY "true")
    endif()
  endif()
endif()
set(FREHG2_GENERATED_DIR "${CMAKE_BINARY_DIR}/generated")
configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/build_info.hpp.in"
  "${FREHG2_GENERATED_DIR}/frehg2/io/build_info.hpp"
  @ONLY)

# ---------------------------------------------------------------------------------------
# Subdirectories
# ---------------------------------------------------------------------------------------
enable_testing()
add_subdirectory(src)
add_subdirectory(tools)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
