# Linear-algebra-layer tests (P2.4 now; P2.5/P2.6 added with the PETSc backend).

# DomainDecomposition tests. These link ONLY frehg2_decomp (Kokkos + MPI, no PETSc).
foreach(t test_decomp_2d test_decomp_3d test_decomp_no_petsc)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::decomp)
  apply_strict_warnings(${t})
endforeach()

# Serial registrations.
add_test(NAME test_decomp_2d COMMAND test_decomp_2d)
add_test(NAME test_decomp_3d COMMAND test_decomp_3d)
add_test(NAME test_decomp_no_petsc COMMAND test_decomp_no_petsc)
set_tests_properties(test_decomp_2d test_decomp_3d test_decomp_no_petsc
  PROPERTIES LABELS "cpu")

# 2-rank registrations for the decompositions.
foreach(t test_decomp_2d test_decomp_3d)
  add_test(NAME ${t}_np2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:${t}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${t}_np2 PROPERTIES LABELS "mpi")
endforeach()

# P2.5/P2.6: PETSc backend, backend-independence, Kokkos->PETSc bridge. These link the
# PETSc backend (PETSc types allowed in tests).
foreach(t test_petsc_linear_solver test_assembly_backend_independent test_kokkos_petsc_bridge)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::petsc_backend)
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
  add_test(NAME ${t}_np2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:${t}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${t}_np2 PROPERTIES LABELS "mpi")
endforeach()
