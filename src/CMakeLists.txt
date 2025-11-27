file(GLOB LIB_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*.c
)
add_executable(${CMAKE_PROJECT_NAME} ${LIB_HEADERS} frehg.cpp)
#   Link OpenMP
IF(NOT Enable_GPU)
    find_package(OpenMP REQUIRED)
    target_link_libraries(frehg PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

#   Find MPI
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     MPI::MPI)

#   Find Kokkos
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE     Kokkos::kokkos)

#   Move the executable into bin/
add_custom_command(TARGET frehg POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:frehg> ../frehg
)
