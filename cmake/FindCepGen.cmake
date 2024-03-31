find_package(PkgConfig)
pkg_check_modules(PC_CepGen QUIET CepGen)

find_path(CepGen_INCLUDE_DIR
    NAMES CepGen/Generator.h
    PATHS ${PC_CepGen_INCLUDE_DIRS} ${CEPGEN_PATH})
find_library(CepGen_LIBRARY
    NAMES CepGen
    PATHS ${PC_CepGEN_LIBRARY_DIRS} ${CEPGEN_PATH}
    PATH_SUFFIXES lib64 lib build)

include(CepGenVersion)
set(CepGen_VERSION ${VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CepGen
    FOUND_VAR CepGen_FOUND
    REQUIRED_VARS CepGen_LIBRARY CepGen_INCLUDE_DIR
    VERSION_VAR CepGen_VERSION)

if(CepGen_FOUND)
    set(CepGen_LIBRARIES ${CepGen_LIBRARY})
    set(CepGen_INCLUDE_DIRS ${CepGen_INCLUDE_DIR})
    set(CepGen_DEFINITIONS ${PC_CepGen_CFLAGS_OTHER})
endif()

if(CepGen_FOUND AND NOT TARGET CepGen::CepGen)
  add_library(CepGen::CepGen UNKNOWN IMPORTED)
  set_target_properties(CepGen::CepGen PROPERTIES
    IMPORTED_LOCATION "${CepGen_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_CepGen_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${CepGen_INCLUDE_DIR}")
endif()

mark_as_advanced(CepGen_INCLUDE_DIR CepGen_LIBRARY)
