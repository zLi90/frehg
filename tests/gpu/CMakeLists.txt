# P10 Task 10.3.5 / 10.3.7 — GPU execution tests.
#
# These are written now but are NEVER run on this macOS machine (no CUDA). They are labeled
# "gpu" and registered DISABLED whenever the build has no CUDA backend (FREHG2_ENABLE_CUDA=OFF),
# so the default local `ctest` skips them and CI never fails because a gpu test cannot execute
# (docs/gpu_validation_policy.md). On a Linux/NVIDIA build configured with -DFREHG2_ENABLE_CUDA=ON
# (and a Kokkos-CUDA + Kokkos-aware PETSc toolchain) they run the real device-path validation.
#
# The bodies are guarded by FREHG2_GPU_ASSEMBLY (GpuAssembly.hpp): on a non-GPU build they
# compile to a trivial skip so the whole suite still builds cleanly here.

# b1-sw end-to-end on GPU (reuses the b1 runner; legacy-reference gate).
add_executable(test_b1_sw_gpu test_b1_sw_gpu.cpp)
target_include_directories(test_b1_sw_gpu PRIVATE ${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_b1_sw_gpu PRIVATE Frehg2::swe Frehg2::petsc_backend)
target_compile_definitions(test_b1_sw_gpu PRIVATE
  FREHG2_LEGACY_B1_DIR="${CMAKE_SOURCE_DIR}/legacy/benchmarks/b1-sw")
apply_strict_warnings(test_b1_sw_gpu)

# b2-gw and async coupling on GPU (Orchestrator path).
foreach(t test_b2_gw_gpu test_async_coupling_gpu)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(${t} PRIVATE Frehg2::orchestrator)
  target_compile_definitions(${t} PRIVATE
    FREHG2_IO_TMP="${CMAKE_CURRENT_BINARY_DIR}"
    FREHG2_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  apply_strict_warnings(${t})
endforeach()

foreach(t test_b1_sw_gpu test_b2_gw_gpu test_async_coupling_gpu)
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "gpu")
  if(NOT FREHG2_ENABLE_CUDA)
    set_tests_properties(${t} PROPERTIES DISABLED TRUE)
  endif()
endforeach()
