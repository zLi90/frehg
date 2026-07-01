# Polygon BC / source-sink tests (P12). Geometry/index/yaml tests are plain host C++; the apply
# test drives real SWE/RE solver state and needs Kokkos (Views). All link Frehg2::bc, which pulls
# in core/io/swe/re transitively. No PETSc is needed (no global solve in these unit tests).
foreach(t test_polygon test_polygon_index test_yaml_polygon test_polygon_apply)
  add_executable(${t} ${t}.cpp)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
  target_link_libraries(${t} PRIVATE Frehg2::bc)
  apply_strict_warnings(${t})
  add_test(NAME ${t} COMMAND ${t})
  set_tests_properties(${t} PROPERTIES LABELS "cpu")
endforeach()
