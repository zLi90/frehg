# P21 performance instrumentation (timing + counters). No PETSc types.
add_library(frehg2_perf STATIC
  Timer.cpp
  PerfRecorder.cpp)
add_library(Frehg2::perf ALIAS frehg2_perf)

target_include_directories(frehg2_perf PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_perf PUBLIC Frehg2::core MPI::MPI_CXX)
apply_strict_warnings(frehg2_perf)
