# FindPETSc.cmake - locate PETSc via its pkg-config file (P1.2.2).
#
# The local PETSc 3.25.1 install ships a pkg-config file (lib/pkgconfig/PETSc.pc) but no
# CMake package config, so this module wraps pkg-config and exposes an imported target
# `PETSc::PETSc`. It also records PETSC_VERSION from petscversion.h.
#
# Hint variables (cache or -D):
#   PETSc_DIR / PETSC_DIR : PETSc install prefix (default /Users/zhili/Codes/local)

if(NOT PETSc_DIR)
  if(DEFINED ENV{PETSC_DIR})
    set(PETSc_DIR "$ENV{PETSC_DIR}")
  elseif(PETSC_DIR)
    set(PETSc_DIR "${PETSC_DIR}")
  else()
    set(PETSc_DIR "/Users/zhili/Codes/local")
  endif()
endif()

find_package(PkgConfig REQUIRED)

# Ensure pkg-config can see PETSc.pc from the install prefix.
set(_petsc_pcdir "${PETSc_DIR}/lib/pkgconfig")
if(EXISTS "${_petsc_pcdir}/PETSc.pc" OR EXISTS "${_petsc_pcdir}/petsc.pc")
  set(ENV{PKG_CONFIG_PATH} "${_petsc_pcdir}:$ENV{PKG_CONFIG_PATH}")
endif()

# PETSc.pc (capital) is the modern name; petsc.pc is the legacy lowercase alias.
pkg_check_modules(PC_PETSC QUIET PETSc)
if(NOT PC_PETSC_FOUND)
  pkg_check_modules(PC_PETSC QUIET petsc)
endif()

# Parse the version from petscversion.h for a robust, pc-independent version string.
set(PETSC_VERSION "")
if(EXISTS "${PETSc_DIR}/include/petscversion.h")
  file(STRINGS "${PETSc_DIR}/include/petscversion.h" _ver_major REGEX "define PETSC_VERSION_MAJOR")
  file(STRINGS "${PETSc_DIR}/include/petscversion.h" _ver_minor REGEX "define PETSC_VERSION_MINOR")
  file(STRINGS "${PETSc_DIR}/include/petscversion.h" _ver_sub REGEX "define PETSC_VERSION_SUBMINOR")
  string(REGEX MATCH "[0-9]+" _vM "${_ver_major}")
  string(REGEX MATCH "[0-9]+" _vm "${_ver_minor}")
  string(REGEX MATCH "[0-9]+" _vs "${_ver_sub}")
  if(_vM AND _vm AND _vs)
    set(PETSC_VERSION "${_vM}.${_vm}.${_vs}")
  endif()
endif()
if(NOT PETSC_VERSION AND PC_PETSC_VERSION)
  set(PETSC_VERSION "${PC_PETSC_VERSION}")
endif()

# Resolve include dir and the petsc library explicitly so we have absolute paths even if
# the .pc only emits -L/-l flags.
find_path(PETSC_INCLUDE_DIR
  NAMES petsc.h
  HINTS ${PC_PETSC_INCLUDE_DIRS} "${PETSc_DIR}/include"
)
find_library(PETSC_LIBRARY
  NAMES petsc
  HINTS ${PC_PETSC_LIBRARY_DIRS} "${PETSc_DIR}/lib"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  REQUIRED_VARS PETSC_LIBRARY PETSC_INCLUDE_DIR
  VERSION_VAR PETSC_VERSION
)

if(PETSc_FOUND AND NOT TARGET PETSc::PETSc)
  add_library(PETSc::PETSc UNKNOWN IMPORTED)
  set_target_properties(PETSc::PETSc PROPERTIES
    IMPORTED_LOCATION "${PETSC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIR}"
  )
  # Carry any extra link flags the .pc reported (e.g. MPI/HDF5 deps) when present.
  if(PC_PETSC_LDFLAGS)
    set_property(TARGET PETSc::PETSc PROPERTY INTERFACE_LINK_LIBRARIES "${PC_PETSC_LDFLAGS}")
  endif()
endif()

mark_as_advanced(PETSC_INCLUDE_DIR PETSC_LIBRARY)
