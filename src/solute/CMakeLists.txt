# Solute / scalar advection-dispersion transport (P8). Standalone solver (wired into the
# Orchestrator in P16). Diffusion goes through the backend-agnostic LinearSolver/SparseSystem
# interface, so there are NO PETSc types in solute/ (enforced by check_solver_seam).
add_library(frehg2_solute STATIC
  SoluteParams.cpp
  Advection.cpp
  Diffusion.cpp
  SourceSink.cpp
  SoluteStepper.cpp)
add_library(Frehg2::solute ALIAS frehg2_solute)

target_include_directories(frehg2_solute PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_solute PUBLIC
  Frehg2::core
  Frehg2::io
  Frehg2::decomp)

apply_strict_warnings(frehg2_solute)
