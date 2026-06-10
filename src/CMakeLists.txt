add_subdirectory(core)
add_subdirectory(io)
add_subdirectory(swe)
add_subdirectory(re)
add_subdirectory(coupling)
add_subdirectory(solute)
add_subdirectory(bc)
add_subdirectory(driver)

add_executable(frehg2 main.cpp)
target_link_libraries(frehg2 PRIVATE frehg2_driver frehg2_solute frehg2_bc)

target_compile_definitions(frehg2 PRIVATE
    $<$<BOOL:${USE_KOKKOS}>:USE_KOKKOS>
    $<$<BOOL:${USE_MPI}>:USE_MPI>
    $<$<BOOL:${USE_PETSC}>:USE_PETSC>
)

if(USE_KOKKOS)
    target_link_libraries(frehg2 PRIVATE Kokkos::kokkos)
endif()

if(USE_MPI)
    target_link_libraries(frehg2 PRIVATE MPI::MPI_CXX)
endif()

if(USE_PETSC)
    target_link_libraries(frehg2 PRIVATE PETSc::PETSc)
endif()
