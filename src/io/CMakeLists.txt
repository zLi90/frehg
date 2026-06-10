add_library(frehg2_io
    AsciiRaster.cpp
    Config.cpp
    Hdf5Reader.cpp
    Hdf5Writer.cpp
    TimeSeries.cpp
)

target_include_directories(frehg2_io
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_io PUBLIC frehg2_core)

target_compile_definitions(frehg2_io
    PUBLIC
        $<$<BOOL:${USE_KOKKOS}>:USE_KOKKOS>
        $<$<BOOL:${USE_MPI}>:USE_MPI>
        $<$<BOOL:${USE_HDF5}>:USE_HDF5>
        $<$<BOOL:${USE_YAML_CPP}>:USE_YAML_CPP>
)

if(USE_HDF5)
    target_link_libraries(frehg2_io PUBLIC ${HDF5_LIBRARIES})
    target_include_directories(frehg2_io PUBLIC ${HDF5_INCLUDE_DIRS})
endif()

if(USE_YAML_CPP)
    target_link_libraries(frehg2_io PUBLIC yaml-cpp::yaml-cpp)
endif()
