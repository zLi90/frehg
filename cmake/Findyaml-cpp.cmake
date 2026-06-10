set(_yamlcpp_roots
    "${YAMLCPP_ROOT}"
    "$ENV{YAMLCPP_ROOT}"
    "/Users/zhili/Codes/local"
)

find_path(yaml-cpp_INCLUDE_DIR
    NAMES yaml-cpp/yaml.h
    PATHS ${_yamlcpp_roots}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

find_library(yaml-cpp_LIBRARY
    NAMES yaml-cpp yaml-cppd
    PATHS ${_yamlcpp_roots}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(yaml-cpp
    REQUIRED_VARS yaml-cpp_INCLUDE_DIR yaml-cpp_LIBRARY
)

if(yaml-cpp_FOUND AND NOT TARGET yaml-cpp::yaml-cpp)
    add_library(yaml-cpp::yaml-cpp UNKNOWN IMPORTED)
    set_target_properties(yaml-cpp::yaml-cpp PROPERTIES
        IMPORTED_LOCATION "${yaml-cpp_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${yaml-cpp_INCLUDE_DIR}"
    )
endif()
