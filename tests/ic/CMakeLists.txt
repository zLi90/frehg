# Flexible initial conditions tests (P14).
add_executable(test_formula_parser test_formula_parser.cpp)
target_include_directories(test_formula_parser PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_formula_parser PRIVATE Frehg2::ic)
apply_strict_warnings(test_formula_parser)
add_test(NAME test_formula_parser COMMAND test_formula_parser)
set_tests_properties(test_formula_parser PROPERTIES LABELS "cpu")

add_executable(test_ic_apply test_ic_apply.cpp)
target_include_directories(test_ic_apply PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_compile_definitions(test_ic_apply PRIVATE
  FREHG2_IO_TMP="${CMAKE_CURRENT_BINARY_DIR}"
  FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_link_libraries(test_ic_apply PRIVATE Frehg2::ic Frehg2::swe Frehg2::re)
apply_strict_warnings(test_ic_apply)
add_test(NAME test_ic_apply COMMAND test_ic_apply)
set_tests_properties(test_ic_apply PROPERTIES LABELS "cpu;regression")
