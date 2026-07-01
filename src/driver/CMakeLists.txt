# Frehg2 production driver (P7). The Orchestrator is the ONLY production path from main() to
# the physics. It is built as a library so the unified driver and the integration tests share
# exactly one code path. It owns no PETSc type: the linear backend is created through the
# backend-agnostic makeLinearSolver() factory (PETSc lives only in frehg2_petsc_backend).
add_library(frehg2_orchestrator STATIC
  Orchestrator.cpp
  AsyncPipeline.cpp)
add_library(Frehg2::orchestrator ALIAS frehg2_orchestrator)

target_include_directories(frehg2_orchestrator PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

# AsyncPipeline (P11) dispatches the groundwater catch-up to a worker thread; link the
# system threads library explicitly so std::thread resolves on the g++-15/libstdc++ toolchain.
find_package(Threads REQUIRED)

target_link_libraries(frehg2_orchestrator PUBLIC
  Frehg2::coupling
  Frehg2::bc
  Frehg2::ic
  Frehg2::monitoring
  Frehg2::solute
  Frehg2::io
  Frehg2::perf
  Frehg2::petsc_backend
  Threads::Threads)

apply_strict_warnings(frehg2_orchestrator)
add_dependencies(frehg2_orchestrator check_no_seqaij check_solver_seam)

add_executable(frehg2 main.cpp)
# orchestrator transitively brings coupling -> swe/re -> core, io, and petsc_backend.
target_link_libraries(frehg2 PRIVATE Frehg2::orchestrator)

apply_strict_warnings(frehg2)

# CQ6 + P2.6.3: building the driver requires the PETSc-pattern and solver-seam checks.
add_dependencies(frehg2 check_no_seqaij check_solver_seam)

set_target_properties(frehg2 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
