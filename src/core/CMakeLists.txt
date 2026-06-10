add_library(frehg2_core
    Domain.cpp
    Grid.cpp
    GwDomain.cpp
    GwState.cpp
    InitialCondition.cpp
    Monitor.cpp
    MpiComm.cpp
    Simulation.cpp
    State.cpp
)

target_include_directories(frehg2_core
    PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/.."
)

target_compile_definitions(frehg2_core
    PUBLIC
        $<$<BOOL:${USE_KOKKOS}>:USE_KOKKOS>
        $<$<BOOL:${USE_MPI}>:USE_MPI>
        $<$<BOOL:${USE_PETSC}>:USE_PETSC>
        $<$<BOOL:${USE_HDF5}>:USE_HDF5>
        $<$<BOOL:${USE_YAML_CPP}>:USE_YAML_CPP>
)

if(USE_KOKKOS)
    target_link_libraries(frehg2_core PUBLIC Kokkos::kokkos)
endif()

if(USE_MPI)
    target_link_libraries(frehg2_core PUBLIC MPI::MPI_CXX)
endif()

if(USE_PETSC)
    target_link_libraries(frehg2_core PUBLIC PETSc::PETSc)
endif()

if(USE_HDF5)
    target_link_libraries(frehg2_core PUBLIC ${HDF5_LIBRARIES})
    target_include_directories(frehg2_core PUBLIC ${HDF5_INCLUDE_DIRS})
endif()

if(USE_YAML_CPP)
    target_link_libraries(frehg2_core PUBLIC yaml-cpp::yaml-cpp)
endif()
