# P21 performance instrumentation tests.

add_executable(test_perf_recorder test_perf_recorder.cpp)
target_include_directories(test_perf_recorder PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(test_perf_recorder PRIVATE Frehg2::perf)
apply_strict_warnings(test_perf_recorder)
add_test(NAME test_perf_recorder COMMAND test_perf_recorder)
set_tests_properties(test_perf_recorder PROPERTIES LABELS "cpu;perf")

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import yaml"
    RESULT_VARIABLE _yaml_rc OUTPUT_QUIET ERROR_QUIET)
  if(_yaml_rc EQUAL 0)
    add_test(NAME test_perf_report_quick
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/perf_report.py
              --quick --bin $<TARGET_FILE:frehg2>
              --out ${CMAKE_CURRENT_BINARY_DIR}/perf_quick)
    set_tests_properties(test_perf_report_quick PROPERTIES
      LABELS "cpu;perf;integration" RESOURCE_LOCK frehg2_run_dirs)
  endif()
endif()
