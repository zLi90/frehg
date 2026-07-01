# P19 unified validation harness tests. Pure-Python checks, registered only when the needed
# interpreters/modules are importable (mirrors tests/io for h5py/PyYAML).

find_package(Python3 COMPONENTS Interpreter QUIET)
if(NOT Python3_FOUND)
  message(STATUS "  validation tests: python3 not found; skipping")
  return()
endif()

# trend_check unit test: pure python, no h5py / frehg2 binary needed.
add_test(NAME test_trend_check
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_trend_check.py)
set_tests_properties(test_trend_check PROPERTIES LABELS "cpu;validation")

# End-to-end harness smoke (capped --quick): runs all 6 benchmarks + b0-lake through the driver
# and asserts the harness completes with no FAIL tier. Requires h5py + PyYAML + the frehg2 binary.
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import h5py, yaml"
  RESULT_VARIABLE _vh_rc OUTPUT_QUIET ERROR_QUIET)
if(_vh_rc EQUAL 0)
  add_test(NAME test_validation_harness_quick
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/run_validation.py
            --quick --bin $<TARGET_FILE:frehg2>
            --out ${CMAKE_CURRENT_BINARY_DIR}/validation_quick)
  # The harness exits non-zero only on a FAIL tier; REVIEW/PASS => 0 (success).
  set_tests_properties(test_validation_harness_quick PROPERTIES
    LABELS "cpu;validation;integration" RESOURCE_LOCK frehg2_run_dirs)
  add_test(NAME test_b0_lake_well_balanced
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/run_b0_lake.py
            --bin $<TARGET_FILE:frehg2> --out ${CMAKE_CURRENT_BINARY_DIR}/b0_lake)
  set_tests_properties(test_b0_lake_well_balanced PROPERTIES
    LABELS "cpu;validation;regression" RESOURCE_LOCK frehg2_run_dirs)
else()
  message(STATUS "  validation tests: h5py/PyYAML not importable; skipping harness smoke")
endif()
