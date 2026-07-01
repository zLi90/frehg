# Polygon-based boundary conditions & source/sink regions (P12). Primary BC path from P12 on;
# the legacy rectangular/bc_type arrays are deprecated. Plain serial/host and free of PETSc
# types — it drives the SWE / RE solvers through their public field accessors (no solver seam
# violation). The plan's separate geometry/ + bcs/ + sources/ dirs are consolidated into this
# single src/bc/ module (they share the PolygonIndex and are tightly coupled).
add_library(frehg2_bc STATIC
  PolygonRegion.cpp
  PolygonIndex.cpp
  PolygonBC.cpp
  PolygonSource.cpp
  PolygonConfig.cpp)
add_library(Frehg2::bc ALIAS frehg2_bc)

target_include_directories(frehg2_bc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_bc PUBLIC
  Frehg2::core
  Frehg2::io
  Frehg2::swe
  Frehg2::re)

apply_strict_warnings(frehg2_bc)
