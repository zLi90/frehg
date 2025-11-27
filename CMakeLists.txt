#   This is the CMakeLists file for building and compiling Frehg

cmake_minimum_required (VERSION 3.20)

# Project Name and Languages
# C is included because legacy parts of frehg likely remain in C
project(frehg LANGUAGES CXX C)

# -----------------------------------------------------------------------------
# 1. Configuration & Standards
# -----------------------------------------------------------------------------
#   If needed, specify the CXX compiler here, otherwise comment out these lines
set(CMAKE_C_COMPILER gcc-14)
set(CMAKE_CXX_COMPILER g++-14)
# Kokkos usually requires at least C++14; C++17 is recommended for modern features.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Enable aggressive warnings to help clean up legacy code
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra>)

# -----------------------------------------------------------------------------
# 2. Dependency: MPI
# -----------------------------------------------------------------------------
find_package(MPI REQUIRED)

# -----------------------------------------------------------------------------
# 3. Dependency: Kokkos
# -----------------------------------------------------------------------------
set(Kokkos_ROOT "/Users/zhili/Codes/kokkos-install")
find_package(Kokkos REQUIRED)

# -----------------------------------------------------------------------------
# 4. Define GPU
# -----------------------------------------------------------------------------
#   Specify if using GPU or not
set(Enable_GPU "OFF")

#   If using GPU, specify the architecture
if (Enable_GPU)
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 80)
    endif()
endif()

#   Add the subdirectory where the source codes are stored
add_subdirectory(src)
