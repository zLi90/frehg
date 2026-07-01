# Frehg2 source tree. Modules are added as they acquire real content per the plan:
#   core   (P1), driver (P1)
#   linear (P2), io (P3), swe (P4), re (P5), coupling (P6), solute (P8),
#   bc (P12), monitoring (P15)
# Placeholder subdirectories carry a CMakeLists.txt now (P1.2.3) and are wired in when
# their sources land.
add_subdirectory(core)
add_subdirectory(linear)
add_subdirectory(io)
add_subdirectory(soil)
add_subdirectory(ic)
add_subdirectory(swe)
add_subdirectory(re)
add_subdirectory(solute)
add_subdirectory(coupling)
add_subdirectory(bc)
add_subdirectory(monitoring)
add_subdirectory(perf)
add_subdirectory(driver)
