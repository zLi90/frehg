# Core data structures (P1.3 types + P2.1/2.2/2.3). Compiled static library.
# Public include roots: include/ (public API headers) and src/ (internal module headers
# included as "core/Grid.hpp", "linear/...").
add_library(frehg2_core STATIC
  Grid.cpp
  Domain.cpp
  GwDomain.cpp
  State.cpp
  GwState.cpp
  MpiComm.cpp)
add_library(Frehg2::core ALIAS frehg2_core)

target_include_directories(frehg2_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_core PUBLIC Kokkos::kokkos MPI::MPI_CXX)

apply_strict_warnings(frehg2_core)
