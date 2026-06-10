add_library(frehg2_coupling
    Coupling.cpp
)

target_include_directories(frehg2_coupling
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_coupling PUBLIC frehg2_core)
