# Solute transport tests (P8). Kernel tests run on Kokkos only; the diffusion/stepper tests
# additionally need a LinearSolver backend (PETSc) for the implicit solve.
foreach(t test_solute_state test_advection test_source_sink test_yaml_solute)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::solute)
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
endforeach()

foreach(t test_diffusion test_solute_stepper)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::solute Frehg2::petsc_backend)
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
endforeach()
