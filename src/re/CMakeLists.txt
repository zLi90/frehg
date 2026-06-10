add_library(frehg2_re
    ReSolver.cpp
)

target_include_directories(frehg2_re
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_link_libraries(frehg2_re PUBLIC frehg2_core frehg2_io)

target_compile_definitions(frehg2_re
    PUBLIC
        $<$<BOOL:${USE_KOKKOS}>:USE_KOKKOS>
        $<$<BOOL:${USE_MPI}>:USE_MPI>
        $<$<BOOL:${USE_PETSC}>:USE_PETSC>
)

if(USE_PETSC)
    target_link_libraries(frehg2_re PUBLIC PETSc::PETSc)
endif()
