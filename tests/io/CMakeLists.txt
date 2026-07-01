# I/O tests (P3). Serial unit tests + MPI multi-rank for the HDF5 io_modes.

# Pure unit tests (no MPI; Config/SHA/TimeSeries need no Kokkos, AsciiRaster inits Kokkos).
foreach(t test_sha256 test_config test_ascii_raster test_timeseries test_yaml_migration)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::io)
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
endforeach()

# test_config needs the benchmark YAML directory.
target_compile_definitions(test_config PRIVATE FREHG2_BENCH_DIR="${CMAKE_SOURCE_DIR}/benchmarks")
# P17 migration test needs the source tree (for legacy v1 fixtures) and benchmark v2 fixtures.
target_compile_definitions(test_yaml_migration PRIVATE
  FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
  FREHG2_BENCH_DIR="${CMAKE_SOURCE_DIR}/benchmarks")

# HDF5 writer/reader/iface/checkpoint tests: MPI + Kokkos; run serial and on 2 ranks.
foreach(t test_hdf5 test_output_writer_iface test_checkpoint)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::io)
  target_compile_definitions(${t} PRIVATE FREHG2_IO_TMP="${CMAKE_CURRENT_BINARY_DIR}")
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
  add_test(NAME ${t}_np2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:${t}> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(${t}_np2 PROPERTIES LABELS "mpi")
endforeach()

# Mark test_hdf5 as producing the h5py-check fixture file (h5py_check.h5).
set_tests_properties(test_hdf5 PROPERTIES FIXTURES_SETUP H5PY_FILE)

# Python/h5py interop check (P3.2 acceptance): only registered if python3 + h5py exist.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_FOUND)
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import h5py"
    RESULT_VARIABLE _h5py_rc OUTPUT_QUIET ERROR_QUIET)
  if(_h5py_rc EQUAL 0)
    add_test(NAME test_hdf5_h5py
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_h5py.py
              ${CMAKE_CURRENT_BINARY_DIR}/h5py_check.h5 4 3)
    set_tests_properties(test_hdf5_h5py PROPERTIES LABELS "cpu" FIXTURES_REQUIRED H5PY_FILE)
  else()
    message(STATUS "  io tests: h5py not importable; skipping test_hdf5_h5py")
  endif()

  # P17: keep docs/yaml_schema_v2.md in sync with the example fixtures (only if PyYAML exists).
  execute_process(COMMAND ${Python3_EXECUTABLE} -c "import yaml"
    RESULT_VARIABLE _yaml_rc OUTPUT_QUIET ERROR_QUIET)
  if(_yaml_rc EQUAL 0)
    add_test(NAME test_yaml_schema_doc_fresh
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/gen_yaml_schema_doc.py --check)
    set_tests_properties(test_yaml_schema_doc_fresh PROPERTIES LABELS "cpu")
  else()
    message(STATUS "  io tests: PyYAML not importable; skipping test_yaml_schema_doc_fresh")
  endif()
else()
  message(STATUS "  io tests: python3 not found; skipping test_hdf5_h5py")
endif()
