# Synchronous SW<->GW coupling tests (P6).

# P6.1 exchange mechanism: sign conventions, K_sat limiting, full-grid coverage. Read-only
# (no solver attach) -> links only the coupling lib.
add_executable(test_coupling_frehg test_coupling_frehg.cpp)
target_include_directories(test_coupling_frehg PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_coupling_frehg PRIVATE Frehg2::coupling)
apply_strict_warnings(test_coupling_frehg)
add_test(NAME test_coupling_frehg COMMAND test_coupling_frehg)
set_tests_properties(test_coupling_frehg PROPERTIES LABELS "cpu")

# P6.1 full-grid per-column independence.
add_executable(test_coupling_full_grid test_coupling_full_grid.cpp)
target_include_directories(test_coupling_full_grid PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_coupling_full_grid PRIVATE Frehg2::coupling)
apply_strict_warnings(test_coupling_full_grid)
add_test(NAME test_coupling_full_grid COMMAND test_coupling_full_grid)
set_tests_properties(test_coupling_full_grid PROPERTIES LABELS "cpu")

# P6.2 synchronous coupled time stepping: exchange conservation + mass-balance/stability.
add_executable(test_coupled_timestep test_coupled_timestep.cpp)
target_include_directories(test_coupled_timestep PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_coupled_timestep PRIVATE Frehg2::coupling Frehg2::petsc_backend)
apply_strict_warnings(test_coupled_timestep)
add_test(NAME test_coupled_timestep COMMAND test_coupled_timestep)
set_tests_properties(test_coupled_timestep PROPERTIES LABELS "cpu;integration")

# P6.3 coupling validation on b1-sw bathymetry with a synthetic soil column.
add_executable(test_coupling_validation test_coupling_validation.cpp)
target_include_directories(test_coupling_validation PRIVATE
  ${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/tests/swe)
target_link_libraries(test_coupling_validation PRIVATE Frehg2::coupling Frehg2::petsc_backend)
target_compile_definitions(test_coupling_validation PRIVATE
  FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
apply_strict_warnings(test_coupling_validation)
add_test(NAME test_coupling_validation COMMAND test_coupling_validation)
set_tests_properties(test_coupling_validation PROPERTIES LABELS "cpu;integration;regression")
