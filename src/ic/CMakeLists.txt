# Flexible initial conditions (P14): formula parser, raster I/O, polygon/constant/formula/restart
# dispatch applied once at Orchestrator::initialize().
add_library(frehg2_ic STATIC
  FormulaParser.cpp
  ICConfig.cpp
  RasterField.cpp
  ICApply.cpp)
add_library(Frehg2::ic ALIAS frehg2_ic)

target_include_directories(frehg2_ic PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_ic PUBLIC
  Frehg2::core
  Frehg2::io
  Frehg2::bc
  Frehg2::swe
  Frehg2::re)

apply_strict_warnings(frehg2_ic)
