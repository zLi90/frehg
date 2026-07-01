# Non-uniform soil tests (P13).
add_executable(test_soil_map test_soil_map.cpp)
target_include_directories(test_soil_map PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_soil_map PRIVATE Frehg2::soil)
apply_strict_warnings(test_soil_map)
add_test(NAME test_soil_map COMMAND test_soil_map)
set_tests_properties(test_soil_map PROPERTIES LABELS "cpu")

add_executable(test_yaml_soil test_yaml_soil.cpp)
target_include_directories(test_yaml_soil PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_yaml_soil PRIVATE Frehg2::soil)
apply_strict_warnings(test_yaml_soil)
add_test(NAME test_yaml_soil COMMAND test_yaml_soil)
set_tests_properties(test_yaml_soil PROPERTIES LABELS "cpu")

# P13 manufactured acceptance gate: per-cell SoilMap dispatch matches per-class uniform runs and
# the classes genuinely differ. Replaces the unusable b4-govindaraju SWE benchmark as the P13
# numerical gate (see test source header).
add_executable(test_nonuniform_soil test_nonuniform_soil.cpp)
target_include_directories(test_nonuniform_soil PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_nonuniform_soil PRIVATE Frehg2::re Frehg2::soil Frehg2::petsc_backend)
apply_strict_warnings(test_nonuniform_soil)
add_test(NAME test_nonuniform_soil COMMAND test_nonuniform_soil)
set_tests_properties(test_nonuniform_soil PROPERTIES LABELS "cpu;regression")
