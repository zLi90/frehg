# Spatially variable soil parameters (P13). A small, Kokkos/PETSc-free data module: a soil-class
# table (SoilParams tuples) + a per-cell class-index raster. The RE solver (Frehg2::re) links this
# and looks up the per-column class; the uniform single-class case reproduces P5 bit-for-bit.
add_library(frehg2_soil STATIC
  SoilMap.cpp
  SoilConfig.cpp)
add_library(Frehg2::soil ALIAS frehg2_soil)

target_include_directories(frehg2_soil PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_soil PUBLIC
  Frehg2::core
  Frehg2::io)

apply_strict_warnings(frehg2_soil)
