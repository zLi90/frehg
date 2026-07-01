# I/O layer (P3): YAML config, ASCII raster, time series, HDF5 writer/reader, checkpoints.
add_library(frehg2_io STATIC
  Sha256.cpp
  Config.cpp
  YamlMigration.cpp
  AsciiRaster.cpp
  TimeSeries.cpp
  Hdf5Support.cpp
  Hdf5Writer.cpp
  Hdf5Reader.cpp
  OutputWriter.cpp)
add_library(Frehg2::io ALIAS frehg2_io)

target_include_directories(frehg2_io PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${FREHG2_GENERATED_DIR}>)

target_include_directories(frehg2_io PUBLIC ${HDF5_C_INCLUDE_DIRS})

target_link_libraries(frehg2_io PUBLIC
  Frehg2::core
  yaml-cpp::yaml-cpp
  ${HDF5_C_LIBRARIES}
  MPI::MPI_CXX)

apply_strict_warnings(frehg2_io)
