# Surface-water solver (P4). Semi-implicit SWE, legacy-exact b1-sw port (stage 4a+).
# Links the backend-agnostic solver layer for matrix assembly (no PETSc types in swe/).
add_library(frehg2_swe STATIC
  SweSolver.cpp)
add_library(Frehg2::swe ALIAS frehg2_swe)

target_include_directories(frehg2_swe PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_swe PUBLIC
  Frehg2::core
  Frehg2::io
  Frehg2::decomp
  Frehg2::perf)

apply_strict_warnings(frehg2_swe)
