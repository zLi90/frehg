cmake_minimum_required(VERSION 3.16)

project(Frehg2 VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(USE_KOKKOS "Enable Kokkos parallel data structures" ON)
option(USE_MPI "Enable MPI support" ON)
option(USE_PETSC "Enable PETSc linear solver support" ON)
option(USE_HDF5 "Enable HDF5 output support" ON)
option(USE_YAML_CPP "Enable yaml-cpp input support" ON)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_DOCS "Build documentation" OFF)

set(FREHG2_LOCAL_PREFIX "/Users/zhili/Codes/local" CACHE PATH "Local dependency prefix")
set(Kokkos_DIR "${FREHG2_LOCAL_PREFIX}/lib/cmake/Kokkos" CACHE PATH "Kokkos CMake package")
set(PETSC_DIR "${FREHG2_LOCAL_PREFIX}" CACHE PATH "PETSc install prefix")
set(MPICH_ROOT "${FREHG2_LOCAL_PREFIX}" CACHE PATH "MPICH install prefix")
set(HDF5_ROOT "${FREHG2_LOCAL_PREFIX}" CACHE PATH "HDF5 install prefix")
set(HDF5_DIR "${FREHG2_LOCAL_PREFIX}/cmake" CACHE PATH "HDF5 CMake package")
set(YAMLCPP_ROOT "${FREHG2_LOCAL_PREFIX}" CACHE PATH "yaml-cpp install prefix")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

if(USE_KOKKOS)
    find_package(Kokkos REQUIRED CONFIG)
    message(STATUS "Kokkos found at ${Kokkos_DIR}")
endif()

if(USE_MPI)
    list(PREPEND CMAKE_PREFIX_PATH "${MPICH_ROOT}")
    find_package(MPI REQUIRED COMPONENTS CXX)
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
endif()

if(USE_PETSC)
    find_package(PETSc REQUIRED)
    message(STATUS "PETSc found at ${PETSC_DIR}")
endif()

if(USE_HDF5)
    find_package(HDF5 COMPONENTS C)
    if(HDF5_FOUND)
        message(STATUS "HDF5 found: ${HDF5_VERSION}")
    else()
        message(WARNING "HDF5 not found under ${HDF5_ROOT}; disabling HDF5 for Phase 1 build.")
        set(USE_HDF5 OFF CACHE BOOL "Enable HDF5 output support" FORCE)
    endif()
endif()

if(USE_YAML_CPP)
    find_package(yaml-cpp MODULE)
    if(yaml-cpp_FOUND)
        message(STATUS "yaml-cpp found")
    else()
        message(WARNING "yaml-cpp not found under ${YAMLCPP_ROOT}; disabling YAML for Phase 1 build.")
        set(USE_YAML_CPP OFF CACHE BOOL "Enable yaml-cpp input support" FORCE)
    endif()
endif()

add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_custom_target(docs
            COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
    endif()
endif()

message(STATUS "")
message(STATUS "Frehg2 Configuration Summary:")
message(STATUS "  CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  CXX Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "  Build Type:   ${CMAKE_BUILD_TYPE}")
message(STATUS "  USE_KOKKOS:   ${USE_KOKKOS}")
message(STATUS "  USE_MPI:      ${USE_MPI}")
message(STATUS "  USE_PETSC:    ${USE_PETSC}")
message(STATUS "  USE_HDF5:     ${USE_HDF5}")
message(STATUS "  USE_YAML_CPP: ${USE_YAML_CPP}")
message(STATUS "  BUILD_TESTS:  ${BUILD_TESTS}")
message(STATUS "")
