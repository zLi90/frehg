# RE solver tests (P5).
add_executable(test_re_constitutive test_re_constitutive.cpp)
target_include_directories(test_re_constitutive PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_re_constitutive PRIVATE Frehg2::re)
apply_strict_warnings(test_re_constitutive)
add_test(NAME test_re_constitutive COMMAND test_re_constitutive)
set_tests_properties(test_re_constitutive PROPERTIES LABELS "cpu")

add_executable(test_re_b2_gw test_re_b2_gw.cpp)
target_include_directories(test_re_b2_gw PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_re_b2_gw PRIVATE Frehg2::re Frehg2::petsc_backend)
target_compile_definitions(test_re_b2_gw PRIVATE
  FREHG2_WARRICK_REF="${CMAKE_SOURCE_DIR}/benchmarks/b2-gw/warrick_water_content_profile.csv"
  FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
apply_strict_warnings(test_re_b2_gw)
add_test(NAME test_re_b2_gw COMMAND test_re_b2_gw)
set_tests_properties(test_re_b2_gw PROPERTIES LABELS "cpu;regression")

# P23: 3D Richards (lateral flux) + per-cell soil — serial gates.
add_executable(test_re_3d test_re_3d.cpp)
target_include_directories(test_re_3d PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_re_3d PRIVATE Frehg2::re Frehg2::soil Frehg2::petsc_backend)
apply_strict_warnings(test_re_3d)
add_test(NAME test_re_3d COMMAND test_re_3d)
set_tests_properties(test_re_3d PROPERTIES LABELS "cpu;regression")

# P23: 3D Richards MPI rank-count equivalence (laterally-coupled). np1 writes the reference.
add_executable(test_re_3d_mpi test_re_3d_mpi.cpp)
target_include_directories(test_re_3d_mpi PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_re_3d_mpi PRIVATE Frehg2::re Frehg2::petsc_backend)
apply_strict_warnings(test_re_3d_mpi)
add_test(NAME test_re_3d_mpi COMMAND test_re_3d_mpi)
set_tests_properties(test_re_3d_mpi PROPERTIES LABELS "cpu;regression" FIXTURES_SETUP re_3d_mpi_ref)
add_test(NAME test_re_3d_mpi_np2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_re_3d_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_re_3d_mpi_np2 PROPERTIES LABELS "mpi;regression" FIXTURES_REQUIRED re_3d_mpi_ref)
add_test(NAME test_re_3d_mpi_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_re_3d_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_re_3d_mpi_np4 PROPERTIES LABELS "mpi;regression" FIXTURES_REQUIRED re_3d_mpi_ref)

# P5.5e: RE MPI rank-count equivalence. np1 writes the reference; np2/np4 compare.
add_executable(test_re_mpi test_re_mpi.cpp)
target_include_directories(test_re_mpi PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_re_mpi PRIVATE Frehg2::re Frehg2::petsc_backend)
apply_strict_warnings(test_re_mpi)
add_test(NAME test_re_mpi COMMAND test_re_mpi)
set_tests_properties(test_re_mpi PROPERTIES LABELS "cpu;regression" FIXTURES_SETUP re_mpi_ref)
add_test(NAME test_re_mpi_np2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_re_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_re_mpi_np2 PROPERTIES LABELS "mpi;regression" FIXTURES_REQUIRED re_mpi_ref)
add_test(NAME test_re_mpi_np4
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_re_mpi> ${MPIEXEC_POSTFLAGS})
set_tests_properties(test_re_mpi_np4 PROPERTIES LABELS "mpi;regression" FIXTURES_REQUIRED re_mpi_ref)
