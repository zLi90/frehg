add_library(frehg2_bc
    BoundaryCondition.cpp
    Polygon.cpp
)

target_include_directories(frehg2_bc
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_bc PUBLIC frehg2_core)
