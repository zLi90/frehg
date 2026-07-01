# Point/line monitor system (P15): probe sampling, line-flux integration, CSV output.
add_library(frehg2_monitoring STATIC
  MonitorConfig.cpp
  ProbeLocator.cpp
  LineFlux.cpp
  MonitorWriter.cpp)
add_library(Frehg2::monitoring ALIAS frehg2_monitoring)

target_include_directories(frehg2_monitoring PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)

target_link_libraries(frehg2_monitoring PUBLIC
  Frehg2::core
  Frehg2::io
  Frehg2::swe
  Frehg2::re)

apply_strict_warnings(frehg2_monitoring)
