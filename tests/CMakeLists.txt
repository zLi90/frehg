add_executable(frehg2_directory_structure_tests
    core/test_directory_structure.cpp
)
target_compile_definitions(frehg2_directory_structure_tests PRIVATE
    FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

add_executable(frehg2_type_tests
    core/test_types.cpp
)
target_link_libraries(frehg2_type_tests PRIVATE frehg2_core)

add_executable(frehg2_grid_tests
    core/test_grid.cpp
)
target_link_libraries(frehg2_grid_tests PRIVATE frehg2_core)

add_executable(frehg2_domain_tests
    core/test_domain.cpp
)
target_link_libraries(frehg2_domain_tests PRIVATE frehg2_core)

add_executable(frehg2_state_tests
    core/test_state.cpp
)
target_link_libraries(frehg2_state_tests PRIVATE frehg2_core)

add_executable(frehg2_ic_tests
    core/test_ic.cpp
)
target_link_libraries(frehg2_ic_tests PRIVATE frehg2_core)

add_executable(frehg2_monitor_tests
    core/test_monitor.cpp
)
target_link_libraries(frehg2_monitor_tests PRIVATE frehg2_core frehg2_bc)

if(USE_MPI)
    add_executable(frehg2_mpi_tests
        core/test_mpi.cpp
    )
    target_link_libraries(frehg2_mpi_tests PRIVATE frehg2_core)
endif()

add_executable(frehg2_config_tests
    io/test_config.cpp
)
target_link_libraries(frehg2_config_tests PRIVATE frehg2_io)
target_compile_definitions(frehg2_config_tests PRIVATE
    FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

add_executable(frehg2_ascii_raster_tests
    io/test_ascii_raster.cpp
)
target_link_libraries(frehg2_ascii_raster_tests PRIVATE frehg2_io)

add_executable(frehg2_timeseries_tests
    io/test_timeseries.cpp
)
target_link_libraries(frehg2_timeseries_tests PRIVATE frehg2_io)

add_executable(frehg2_hdf5_tests
    io/test_hdf5.cpp
)
target_link_libraries(frehg2_hdf5_tests PRIVATE frehg2_io)

add_executable(frehg2_1d_swe_matrix_tests
    swe/test_1d_swe_matrix.cpp
)
target_link_libraries(frehg2_1d_swe_matrix_tests PRIVATE frehg2_swe)

add_executable(frehg2_swe_init_tests
    swe/test_swe_init.cpp
)
target_link_libraries(frehg2_swe_init_tests PRIVATE frehg2_swe)

add_executable(frehg2_swe_solver_tests
    swe/test_swe_solver.cpp
)
target_link_libraries(frehg2_swe_solver_tests PRIVATE frehg2_swe)

add_executable(frehg2_swe_sourcesink_tests
    swe/test_swe_sourcesink.cpp
)
target_link_libraries(frehg2_swe_sourcesink_tests PRIVATE frehg2_swe)

add_executable(frehg2_swe_cfl_tests
    swe/test_swe_cfl.cpp
)
target_link_libraries(frehg2_swe_cfl_tests PRIVATE frehg2_swe)

add_executable(frehg2_1d_re_column_tests
    re/test_1d_re_column.cpp
)
target_link_libraries(frehg2_1d_re_column_tests PRIVATE frehg2_re)

add_executable(frehg2_re_init_tests
    re/test_re_init.cpp
)
target_link_libraries(frehg2_re_init_tests PRIVATE frehg2_re)

add_executable(frehg2_re_matrix_tests
    re/test_re_matrix.cpp
)
target_link_libraries(frehg2_re_matrix_tests PRIVATE frehg2_re)

add_executable(frehg2_re_solver_tests
    re/test_re_solver.cpp
)
target_link_libraries(frehg2_re_solver_tests PRIVATE frehg2_re)

add_executable(frehg2_re_sourcesink_tests
    re/test_re_sourcesink.cpp
)
target_link_libraries(frehg2_re_sourcesink_tests PRIVATE frehg2_re)

add_executable(frehg2_nonuniform_soil_tests
    re/test_nonuniform_soil.cpp
)
target_link_libraries(frehg2_nonuniform_soil_tests PRIVATE frehg2_re)

add_executable(frehg2_soil_map_tests
    re/test_soil_map.cpp
)
target_link_libraries(frehg2_soil_map_tests PRIVATE frehg2_re)
target_compile_definitions(frehg2_soil_map_tests PRIVATE
    FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

add_executable(frehg2_coupling_frehg_tests
    coupling/test_coupling_frehg.cpp
)
target_link_libraries(frehg2_coupling_frehg_tests PRIVATE frehg2_coupling)

add_executable(frehg2_coupled_timestep_tests
    coupling/test_coupled_timestep.cpp
)
target_link_libraries(frehg2_coupled_timestep_tests PRIVATE frehg2_coupling)

add_executable(frehg2_async_coupling_tests
    coupling/test_async_coupling.cpp
)
target_link_libraries(frehg2_async_coupling_tests PRIVATE frehg2_coupling)

add_executable(frehg2_different_dt_tests
    coupling/test_different_dt.cpp
)
target_link_libraries(frehg2_different_dt_tests PRIVATE frehg2_coupling)

add_executable(frehg2_polygon_tests
    bc/test_polygon.cpp
)
target_link_libraries(frehg2_polygon_tests PRIVATE frehg2_bc)

add_executable(frehg2_advanced_bc_tests
    bc/test_advanced_bc.cpp
)
target_link_libraries(frehg2_advanced_bc_tests PRIVATE frehg2_bc)

add_executable(frehg2_polygon_sourcesink_tests
    bc/test_sourcesink.cpp
)
target_link_libraries(frehg2_polygon_sourcesink_tests PRIVATE frehg2_bc)

add_executable(frehg2_benchmark_polygon_config_tests
    bc/test_benchmark_polygon_config.cpp
)
target_link_libraries(frehg2_benchmark_polygon_config_tests PRIVATE frehg2_bc frehg2_io)
target_compile_definitions(frehg2_benchmark_polygon_config_tests PRIVATE
    FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

add_executable(frehg2_solute_solver_tests
    solute/test_solute_solver.cpp
)
target_link_libraries(frehg2_solute_solver_tests PRIVATE frehg2_solute)

add_executable(frehg2_solute_coupling_tests
    solute/test_solute_coupling.cpp
)
target_link_libraries(frehg2_solute_coupling_tests PRIVATE frehg2_solute)

add_test(NAME frehg2_directory_structure_tests COMMAND frehg2_directory_structure_tests)
add_test(NAME frehg2_type_tests COMMAND frehg2_type_tests)
add_test(NAME frehg2_grid_tests COMMAND frehg2_grid_tests)
add_test(NAME frehg2_domain_tests COMMAND frehg2_domain_tests)
add_test(NAME frehg2_state_tests COMMAND frehg2_state_tests)
add_test(NAME frehg2_ic_tests COMMAND frehg2_ic_tests)
add_test(NAME frehg2_monitor_tests COMMAND frehg2_monitor_tests)
if(USE_MPI)
    add_test(
        NAME frehg2_mpi_tests
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
                "$<TARGET_FILE:frehg2_mpi_tests>"
    )
endif()
add_test(NAME frehg2_config_tests COMMAND frehg2_config_tests)
add_test(NAME frehg2_ascii_raster_tests COMMAND frehg2_ascii_raster_tests)
add_test(NAME frehg2_timeseries_tests COMMAND frehg2_timeseries_tests)
add_test(NAME frehg2_hdf5_tests COMMAND frehg2_hdf5_tests)
add_test(NAME frehg2_1d_swe_matrix_tests COMMAND frehg2_1d_swe_matrix_tests)
add_test(NAME frehg2_swe_init_tests COMMAND frehg2_swe_init_tests)
add_test(NAME frehg2_swe_solver_tests COMMAND frehg2_swe_solver_tests)
add_test(NAME frehg2_swe_sourcesink_tests COMMAND frehg2_swe_sourcesink_tests)
add_test(NAME frehg2_swe_cfl_tests COMMAND frehg2_swe_cfl_tests)
add_test(NAME frehg2_1d_re_column_tests COMMAND frehg2_1d_re_column_tests)
add_test(NAME frehg2_re_init_tests COMMAND frehg2_re_init_tests)
add_test(NAME frehg2_re_matrix_tests COMMAND frehg2_re_matrix_tests)
add_test(NAME frehg2_re_solver_tests COMMAND frehg2_re_solver_tests)
add_test(NAME frehg2_re_sourcesink_tests COMMAND frehg2_re_sourcesink_tests)
add_test(NAME frehg2_nonuniform_soil_tests COMMAND frehg2_nonuniform_soil_tests)
add_test(NAME frehg2_soil_map_tests COMMAND frehg2_soil_map_tests)
add_test(NAME frehg2_coupling_frehg_tests COMMAND frehg2_coupling_frehg_tests)
add_test(NAME frehg2_coupled_timestep_tests COMMAND frehg2_coupled_timestep_tests)
add_test(NAME frehg2_async_coupling_tests COMMAND frehg2_async_coupling_tests)
add_test(NAME frehg2_different_dt_tests COMMAND frehg2_different_dt_tests)
add_test(NAME frehg2_polygon_tests COMMAND frehg2_polygon_tests)
add_test(NAME frehg2_advanced_bc_tests COMMAND frehg2_advanced_bc_tests)
add_test(NAME frehg2_polygon_sourcesink_tests COMMAND frehg2_polygon_sourcesink_tests)
add_test(NAME frehg2_benchmark_polygon_config_tests COMMAND frehg2_benchmark_polygon_config_tests)
add_test(NAME frehg2_solute_solver_tests COMMAND frehg2_solute_solver_tests)
add_test(NAME frehg2_solute_coupling_tests COMMAND frehg2_solute_coupling_tests)
add_test(NAME frehg2_help COMMAND frehg2 --help)
add_test(NAME frehg2_placeholder_run COMMAND frehg2 "${CMAKE_SOURCE_DIR}/benchmarks/b1-sw/b1-sw.yaml")
add_test(
    NAME frehg2_b1_production_schema_run
    COMMAND frehg2 "${CMAKE_SOURCE_DIR}/benchmarks/b1-sw/b1-sw.production.yaml"
)
add_test(
    NAME frehg2_b2_production_schema_run
    COMMAND frehg2 "${CMAKE_SOURCE_DIR}/benchmarks/b2-gw/b2-gw.production.yaml"
)
add_test(
    NAME frehg2_coupled_smoke_run
    COMMAND frehg2 "${CMAKE_SOURCE_DIR}/benchmarks/p15-coupled-smoke/p15-coupled-smoke.yaml"
)
add_test(
    NAME frehg2_p16_runtime_surface_run
    COMMAND frehg2 "${CMAKE_SOURCE_DIR}/benchmarks/p16-runtime-surface/p16-runtime-surface.yaml"
)
add_test(
    NAME frehg2_p16_runtime_groundwater_run
    COMMAND frehg2 "${CMAKE_SOURCE_DIR}/benchmarks/p16-runtime-groundwater/p16-runtime-groundwater.yaml"
)
add_test(
    NAME frehg2_serghei_converter_b3
    COMMAND python3 "${CMAKE_SOURCE_DIR}/scripts/serghei_to_frehg2.py"
            --input "${CMAKE_SOURCE_DIR}/legacy/benchmarks/b3-kirkland"
            --output "${CMAKE_BINARY_DIR}/p17-conversion/b3-kirkland"
            --case-id b3-kirkland
            --force
)
add_test(
    NAME frehg2_serghei_converter_b4
    COMMAND python3 "${CMAKE_SOURCE_DIR}/scripts/serghei_to_frehg2.py"
            --input "${CMAKE_SOURCE_DIR}/legacy/benchmarks/b4-govindaraju"
            --output "${CMAKE_BINARY_DIR}/p17-conversion/b4-govindaraju"
            --case-id b4-govindaraju
            --force
)
add_test(
    NAME frehg2_serghei_converter_b5
    COMMAND python3 "${CMAKE_SOURCE_DIR}/scripts/serghei_to_frehg2.py"
            --input "${CMAKE_SOURCE_DIR}/legacy/benchmarks/b5-vcatchment"
            --output "${CMAKE_BINARY_DIR}/p17-conversion/b5-vcatchment"
            --case-id b5-vcatchment
            --force
)
add_test(
    NAME frehg2_serghei_converter_b6
    COMMAND python3 "${CMAKE_SOURCE_DIR}/scripts/serghei_to_frehg2.py"
            --input "${CMAKE_SOURCE_DIR}/legacy/benchmarks/b6-superslab"
            --output "${CMAKE_BINARY_DIR}/p17-conversion/b6-superslab"
            --case-id b6-superslab
            --force
)
add_test(
    NAME frehg2_validation_runner_converted_cases
    COMMAND python3 "${CMAKE_SOURCE_DIR}/scripts/run_benchmarks.py"
            --repo-root "${CMAKE_SOURCE_DIR}"
            --cases b3-kirkland b4-govindaraju b5-vcatchment b6-superslab
            --report-dir "${CMAKE_BINARY_DIR}/p18-validation-converted"
)

set_tests_properties(
    frehg2_domain_tests
    frehg2_state_tests
    frehg2_ic_tests
    frehg2_monitor_tests
    frehg2_swe_init_tests
    frehg2_1d_re_column_tests
    frehg2_re_init_tests
    frehg2_re_matrix_tests
    frehg2_re_solver_tests
    frehg2_re_sourcesink_tests
    frehg2_nonuniform_soil_tests
    frehg2_soil_map_tests
    frehg2_coupling_frehg_tests
    frehg2_coupled_timestep_tests
    frehg2_async_coupling_tests
    frehg2_different_dt_tests
    frehg2_polygon_tests
    frehg2_advanced_bc_tests
    frehg2_polygon_sourcesink_tests
    frehg2_benchmark_polygon_config_tests
    frehg2_solute_solver_tests
    frehg2_solute_coupling_tests
    frehg2_placeholder_run
    frehg2_b1_production_schema_run
    frehg2_b2_production_schema_run
    frehg2_coupled_smoke_run
    frehg2_p16_runtime_surface_run
    frehg2_p16_runtime_groundwater_run
    frehg2_serghei_converter_b3
    frehg2_serghei_converter_b4
    frehg2_serghei_converter_b5
    frehg2_serghei_converter_b6
    frehg2_validation_runner_converted_cases
    PROPERTIES ENVIRONMENT "OMP_PROC_BIND=false"
)
