add_library(frehg2_swe
    SweSolver.cpp
)

target_include_directories(frehg2_swe
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_swe PUBLIC frehg2_core frehg2_io)

target_compile_definitions(frehg2_swe
    PUBLIC
        $<$<BOOL:${USE_KOKKOS}>:USE_KOKKOS>
        $<$<BOOL:${USE_MPI}>:USE_MPI>
        $<$<BOOL:${USE_PETSC}>:USE_PETSC>
)

if(USE_PETSC)
    target_link_libraries(frehg2_swe PUBLIC PETSc::PETSc)
endif()
