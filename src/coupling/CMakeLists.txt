# Synchronous SW<->GW coupling (P6). Drives SweSolver/ReSolver through their public state;
# no PETSc types here (the solver/linear seam stays inside each solver). Async coupling P11.
add_library(frehg2_coupling STATIC
  Coupling.cpp)
add_library(Frehg2::coupling ALIAS frehg2_coupling)

target_include_directories(frehg2_coupling PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_coupling PUBLIC
  Frehg2::core
  Frehg2::swe
  Frehg2::re)

apply_strict_warnings(frehg2_coupling)
