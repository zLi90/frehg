add_library(frehg2_driver
    SimulationDriver.cpp
)

target_include_directories(frehg2_driver
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_driver
    PUBLIC
        frehg2_coupling
        frehg2_io
        frehg2_re
        frehg2_swe
)
