# Developer tools (not part of the production executable). P2B.2 solver bake-off harness.
add_executable(solver_bakeoff solver_bakeoff.cpp)
target_link_libraries(solver_bakeoff PRIVATE Frehg2::petsc_backend)
apply_strict_warnings(solver_bakeoff)
set_target_properties(solver_bakeoff PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

add_executable(gen_b1_np1_reference gen_b1_np1_reference.cpp)
target_include_directories(gen_b1_np1_reference PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(gen_b1_np1_reference PRIVATE Frehg2::swe Frehg2::petsc_backend)
target_compile_definitions(gen_b1_np1_reference PRIVATE
  FREHG2_LEGACY_B1_DIR="${CMAKE_SOURCE_DIR}/legacy/benchmarks/b1-sw")
apply_strict_warnings(gen_b1_np1_reference)
set_target_properties(gen_b1_np1_reference PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

add_executable(gen_b1_reference gen_b1_reference.cpp)
target_include_directories(gen_b1_reference PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(gen_b1_reference PRIVATE Frehg2::swe Frehg2::petsc_backend)
apply_strict_warnings(gen_b1_reference)
set_target_properties(gen_b1_reference PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

# P17: v1/experimental -> frozen v2 YAML migration tool.
add_executable(migrate_yaml_v1_to_v2 migrate_yaml_v1_to_v2.cpp)
target_link_libraries(migrate_yaml_v1_to_v2 PRIVATE Frehg2::io)
apply_strict_warnings(migrate_yaml_v1_to_v2)
set_target_properties(migrate_yaml_v1_to_v2 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
