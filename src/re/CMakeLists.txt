# Richards-equation solver (P5). PCA head solve + flux corrector; no PETSc in re/.
add_library(frehg2_re STATIC
  ReSolver.cpp)
add_library(Frehg2::re ALIAS frehg2_re)

target_include_directories(frehg2_re PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_re PUBLIC
  Frehg2::core
  Frehg2::io
  Frehg2::soil
  Frehg2::decomp
  Frehg2::perf)

apply_strict_warnings(frehg2_re)
