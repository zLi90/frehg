# Core tests (P1.3, P1.5).

# P1.3: core type definitions.
add_executable(test_types test_types.cpp)
target_include_directories(test_types PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_types PRIVATE Frehg2::core)
apply_strict_warnings(test_types)
add_test(NAME test_types COMMAND test_types)
set_tests_properties(test_types PROPERTIES LABELS "cpu")

# P2.1/2.2: Grid, Domain/GwDomain, State/GwState (header + core lib).
foreach(t test_grid test_domain test_state)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::core)
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
endforeach()

# P9: on-node parallel-loop helpers. Verifies parallelForRange/Surface/Volume reproduce serial
# element-wise writes and that Kokkos Max/Min reductions match the serial scan bit-for-bit.
add_executable(test_parallel_for test_parallel_for.cpp)
target_include_directories(test_parallel_for PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_parallel_for PRIVATE Frehg2::core)
apply_strict_warnings(test_parallel_for)
add_test(NAME test_parallel_for COMMAND test_parallel_for)
set_tests_properties(test_parallel_for PROPERTIES LABELS "cpu")
# GPU host/device equivalence stub: same binary, labeled "gpu" and DISABLED (no CUDA locally).
# On a Kokkos-CUDA build, remove DISABLED to run the wrappers on the device (P10 validation).
add_test(NAME test_parallel_for_gpu COMMAND test_parallel_for)
set_tests_properties(test_parallel_for_gpu PROPERTIES LABELS "gpu" DISABLED TRUE)

# P2.3: MpiComm. Runs serial (direct) and on 4 ranks (local MPICH mpiexec).
add_executable(test_mpi test_mpi.cpp)
target_include_directories(test_mpi PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_mpi PRIVATE Frehg2::core)
apply_strict_warnings(test_mpi)
add_test(NAME test_mpi COMMAND test_mpi)
set_tests_properties(test_mpi PROPERTIES LABELS "cpu")
add_test(NAME test_mpi_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_mpi_np4 PROPERTIES LABELS "mpi")

# P1.6: verify the local PETSc supports the MATMPIAIJ solver path (and report default).
add_executable(test_petsc_mat test_petsc_mat.cpp)
target_link_libraries(test_petsc_mat PRIVATE MPI::MPI_CXX PETSc::PETSc)
apply_strict_warnings(test_petsc_mat)
add_test(NAME test_petsc_mat COMMAND test_petsc_mat)
set_tests_properties(test_petsc_mat PROPERTIES LABELS "cpu")

# CQ1: strict warnings must REJECT non-warning-clean code. Compile the bad fixture with
# the exact strict flags; WILL_FAIL means the test passes iff compilation fails.
add_test(NAME cq1_strict_warnings_rejects_bad_code
  COMMAND ${CMAKE_CXX_COMPILER}
          -std=c++20 -Wall -Wextra -Werror -Wno-unused-parameter
          -c ${CMAKE_CURRENT_SOURCE_DIR}/test_cq1_strict_warnings.cpp
          -o ${CMAKE_CURRENT_BINARY_DIR}/cq1_bad.o)
set_tests_properties(cq1_strict_warnings_rejects_bad_code PROPERTIES
  WILL_FAIL TRUE
  LABELS "cpu")
