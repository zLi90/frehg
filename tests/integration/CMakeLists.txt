# P7 multi-module integration tests: the unified Orchestrator production path.
# Each test links the orchestrator library (which transitively pulls coupling/swe/re/io and
# the PETSc backend) and uses an in-binary-dir scratch path for output/checkpoints.

function(frehg2_add_integration_test name)
  add_executable(${name} ${name}.cpp)
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${name} PRIVATE Frehg2::orchestrator)
  target_compile_definitions(${name} PRIVATE
    FREHG2_IO_TMP="${CMAKE_CURRENT_BINARY_DIR}"
    FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  apply_strict_warnings(${name})
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS "cpu;integration")
endfunction()

frehg2_add_integration_test(test_smoke_all_modes)
frehg2_add_integration_test(test_grid_from_yaml)

# P17 production schema gate + v1->v2 migration round-trip on the driver.
frehg2_add_integration_test(test_schema_gate)
# Capability gate: fail fast on parsed-but-unimplemented options (wind/subgrid/diffusive_wave/
# baroclinic, evap models, lateral GW BCs, iterative GW solves, bad output.variables).
frehg2_add_integration_test(test_capability_gate)
frehg2_add_integration_test(test_coupled_timestepping)
frehg2_add_integration_test(test_simulation_summary)

# Parity + restart compare against the direct solver paths -> also regression-labelled.
frehg2_add_integration_test(test_orchestrator_parity)
set_tests_properties(test_orchestrator_parity PROPERTIES LABELS "cpu;integration;regression")

# P11 async-vs-sync coupling equivalence (the async gate; bit-identical to synchronous).
frehg2_add_integration_test(test_async_coupling)
set_tests_properties(test_async_coupling PROPERTIES LABELS "cpu;integration;regression")
frehg2_add_integration_test(test_restart)
set_tests_properties(test_restart PROPERTIES LABELS "cpu;integration;regression")

# P12 polygon BC / source-sink gate: channel outflow matches the prescribed discharge.
frehg2_add_integration_test(test_polygon_bc)
set_tests_properties(test_polygon_bc PROPERTIES LABELS "cpu;integration;regression")

# P13 non-uniform soil through the Orchestrator: soil.map is parsed/sliced/applied via run().
frehg2_add_integration_test(test_soil_orchestrator)
set_tests_properties(test_soil_orchestrator PROPERTIES LABELS "cpu;integration;regression")

# P14 flexible IC: YAML restart type matches the CLI --restart path.
frehg2_add_integration_test(test_ic_restart_yaml)
set_tests_properties(test_ic_restart_yaml PROPERTIES LABELS "cpu;integration;regression")

frehg2_add_integration_test(test_monitors_orchestrator)
set_tests_properties(test_monitors_orchestrator PROPERTIES LABELS "cpu;integration;regression")

# P16 solute production integration: stepper wired into Orchestrator::step(); mass balance,
# monitor "C" column, restart preserves conc.
frehg2_add_integration_test(test_solute_in_orchestrator)
set_tests_properties(test_solute_in_orchestrator PROPERTIES LABELS "cpu;integration;regression")
frehg2_add_integration_test(test_solute_restart)
set_tests_properties(test_solute_restart PROPERTIES LABELS "cpu;integration;regression")

# P16-completion: SW<->GW solute exchange wired to the coupling water exchange. A coupled,
# quiescent, closed column transports a surface solute slug into the groundwater and conserves
# total solute mass to machine precision.
frehg2_add_integration_test(test_solute_exchange)
set_tests_properties(test_solute_exchange PROPERTIES LABELS "cpu;integration;regression")

# P18 b4-b6 SERGHEI/legacy benchmark ports (review-tier; capped-step runs through the driver).
foreach(b test_b4_govindaraju test_b5_vcatchment test_b6_kuan)
  frehg2_add_integration_test(${b})
  set_tests_properties(${b} PROPERTIES LABELS "cpu;integration;benchmark")
endforeach()

# b3-kirkland: FULL 2-D layered infiltration (P23 3-D Richards + per-cell soil + fixed-flux top).
# Runs to t_end (~1 min, ~8.7k adaptive steps) and compares the h=0 / h=-400 head contours to the
# digitized Kirkland reference (review-tier). Longer timeout than the capped b4-b6 ports.
frehg2_add_integration_test(test_b3_kirkland)
set_tests_properties(test_b3_kirkland PROPERTIES LABELS "cpu;integration;benchmark" TIMEOUT 600)
