# Surface-water solver tests (P4). Stage 4a: initialization + face geometry.
# Serial CPU unit tests (SWE solver is plain serial/host through Phase 4; DP3 defers
# on-node Kokkos/GPU to P9/P10).
foreach(t test_swe_init test_swe_geometry)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::swe)
  target_compile_definitions(${t} PRIVATE FREHG2_BENCH_DIR="${CMAKE_SOURCE_DIR}/benchmarks")
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
endforeach()

# P4.0 + 4c: benchmark gates. These need the concrete one-rank PETSc solver behind the
# LinearSolver interface (PETSc types allowed in tests).
add_executable(test_swe_b0 test_swe_b0.cpp)
target_include_directories(test_swe_b0 PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_swe_b0 PRIVATE Frehg2::swe Frehg2::petsc_backend)
apply_strict_warnings(test_swe_b0)
add_test(NAME test_swe_b0 COMMAND test_swe_b0)
set_tests_properties(test_swe_b0 PROPERTIES LABELS "cpu;regression")

add_executable(test_swe_b1 test_swe_b1.cpp)
target_include_directories(test_swe_b1 PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_swe_b1 PRIVATE Frehg2::swe Frehg2::petsc_backend)
target_compile_definitions(test_swe_b1 PRIVATE
  FREHG2_LEGACY_B1_DIR="${CMAKE_SOURCE_DIR}/legacy/benchmarks/b1-sw")
apply_strict_warnings(test_swe_b1)
add_test(NAME test_swe_b1 COMMAND test_swe_b1)
set_tests_properties(test_swe_b1 PROPERTIES LABELS "cpu;regression")

add_executable(test_swe_b1_mpi test_swe_b1_mpi.cpp)
target_include_directories(test_swe_b1_mpi PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_swe_b1_mpi PRIVATE Frehg2::swe Frehg2::petsc_backend)
target_compile_definitions(test_swe_b1_mpi PRIVATE
  FREHG2_LEGACY_B1_DIR="${CMAKE_SOURCE_DIR}/legacy/benchmarks/b1-sw"
  FREHG2_B1_NP1_REF="${CMAKE_SOURCE_DIR}/benchmarks/b1-sw/frehg2_np1_depth_reference.txt")
apply_strict_warnings(test_swe_b1_mpi)
add_test(NAME test_swe_b1_mpi_np2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_swe_b1_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_swe_b1_mpi_np2 PROPERTIES LABELS "mpi;regression")
add_test(NAME test_swe_b1_mpi_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_swe_b1_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_swe_b1_mpi_np4 PROPERTIES LABELS "mpi;regression")
