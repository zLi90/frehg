add_library(frehg2_solute
    SoluteSolver.cpp
)

target_include_directories(frehg2_solute
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_solute PUBLIC frehg2_core)
