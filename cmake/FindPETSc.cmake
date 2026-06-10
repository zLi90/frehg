find_package(PkgConfig QUIET)

set(_petsc_roots
    "${PETSC_DIR}"
    "$ENV{PETSC_DIR}"
    "/Users/zhili/Codes/local"
)

if(PkgConfig_FOUND)
    foreach(_root IN LISTS _petsc_roots)
        if(_root)
            list(APPEND _petsc_pkg_paths "${_root}/lib/pkgconfig")
        endif()
    endforeach()
    list(JOIN _petsc_pkg_paths ":" _petsc_pkg_path)
    set(ENV{PKG_CONFIG_PATH} "${_petsc_pkg_path}:$ENV{PKG_CONFIG_PATH}")
    pkg_check_modules(PC_PETSC QUIET PETSc)
endif()

find_path(PETSC_INCLUDE_DIR
    NAMES petsc.h
    HINTS ${PC_PETSC_INCLUDE_DIRS}
    PATHS ${_petsc_roots}
    PATH_SUFFIXES include
)

find_library(PETSC_LIBRARY
    NAMES petsc
    HINTS ${PC_PETSC_LIBRARY_DIRS}
    PATHS ${_petsc_roots}
    PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSC_INCLUDE_DIR PETSC_LIBRARY
)

if(PETSc_FOUND AND NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc UNKNOWN IMPORTED)
    set_target_properties(PETSc::PETSc PROPERTIES
        IMPORTED_LOCATION "${PETSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIR}"
    )
endif()

set(PETSc_INCLUDE_DIRS "${PETSC_INCLUDE_DIR}")
set(PETSc_LIBRARIES "${PETSC_LIBRARY}")
