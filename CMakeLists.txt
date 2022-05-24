#--- first check if we are at CERN
set(LXPLUS_ENV "${CMAKE_SOURCE_DIR}/source-lxplus.sh")
if($ENV{HOSTNAME} MATCHES "^lxplus[0-9]+.cern.ch")
  if(NOT DEFINED ENV{CEPGEN_LXPLUS_ENV})
    message(FATAL_ERROR "LXPLUS environment variables not set!\nPlease source ${LXPLUS_ENV}\nand clean up this build environment, i.e.\nrm -rf CMake*\nand try again...")
  endif()
  message(STATUS "Compiling on LXPLUS.")
  set(IS_LXPLUS "yes")
endif()

set(CPACK_RPM_COMPONENT_INSTALL ON)
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)
set(VERSION 1.0.0)

option(CMAKE_BUILD_CORE "Build Core library" ON)
option(CMAKE_BUILD_PROCESSES "Build physics processes" ON)
option(CMAKE_BUILD_ADDONS "Build AddOns" ON)
option(CMAKE_BUILD_SCAN "Build cross-section scan utilitary" OFF)
option(CMAKE_BUILD_CARD_CONV "Build steering cards converter utilitary" OFF)
option(CMAKE_BUILD_TESTS "Build tests" OFF)
option(CMAKE_BUILD_FOAM "Build FOAM integrator" OFF)
option(CMAKE_BUILD_UTILS "Build miscellaneous utilities" ON)
option(CMAKE_COVERAGE "Generate code coverage" OFF)
#----- release build by default

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#----- check that at least one target is to be built

if(NOT CMAKE_BUILD_CORE AND NOT CMAKE_BUILD_PROCESSES AND NOT CMAKE_BUILD_ADDONS)
  message(FATAL_ERROR "Nothing to be built!")
endif()

if(CMAKE_BUILD_TESTS)
  enable_testing()
  include(CTest)
  if(CMAKE_COVERAGE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(LINK_FLAGS "${LINK_FLAGS} -fprofile-arcs -ftest-coverage")
    add_custom_command(OUTPUT run_coverage
      COMMAND ctest
      COMMAND lcov --capture --directory CepGen --output-file coverage.info
      COMMAND genhtml -q coverage.info --output-directory coverage
      COMMENT Collecting and creating coverage information)
    add_custom_target(coverage DEPENDS run_coverage)
  endif()
endif()

#set(CMAKE_INSTALL_RPATH "$ORIGIN")
#set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "SET DEFAULT INSTALL_PREFIX" FORCE)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake $ENV{ROOTSYS}/cmake)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(CMAKE_CXX_FLAGS_DEBUG "-pg") # for gprof
#set(CMAKE_CXX_FLAGS_DEBUG "-fprofile-arcs -ftest-coverage") # for gcov
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -O2")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(ROOT_CXX_STANDARD 17)

#----- include external paths

include(GNUInstallDirs)
include(UseEnvironment)

set(CEPGEN_LIBRARIES)
set(CEPGEN_LIB_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_CORE_EXT ${GSL_LIB} ${GSL_CBLAS_LIB} dl)
set(CEPGEN_HAS_FUNCTIONAL OFF)
set(CEPGEN_ADDONS_FILE "${CMAKE_CURRENT_BINARY_DIR}/CepGenAddOns.txt")

include(BuildUtils)

#----- enable fortran for external libraries linking

enable_language(Fortran)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
if(CMAKE_BUILD_PROCESSES)
  add_subdirectory(CepGenProcesses)
endif()
if(CMAKE_BUILD_ADDONS)
  file(REMOVE ${CEPGEN_ADDONS_FILE})
  add_subdirectory(CepGenAddOns)
endif()
if(CMAKE_BUILD_CORE)
  add_subdirectory(CepGen)
  list(APPEND CEPGEN_LIBRARIES CepGen)
endif()

#----- copy the input cards and other files

file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} Cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
file(GLOB external_files ${PROJECT_SOURCE_DIR}/External/*.dat ${PROJECT_SOURCE_DIR}/External/*.mcd)
file(COPY ${external_files} DESTINATION ${PROJECT_BINARY_DIR})
file(GLOB readme_file ${PROJECT_SOURCE_DIR}/CepGen/README)
file(COPY ${readme_file} DESTINATION ${PROJECT_BINARY_DIR})

#----- set the tests/utils directory

if(CMAKE_BUILD_CORE)
  #--- build the executables
  add_subdirectory(src)
  #--- build the tests
  if(CMAKE_BUILD_TESTS)
    add_subdirectory(test)
  endif()
endif()

#----- installation rules

install(DIRECTORY CepGen
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT devel FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Cards
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/CepGen
  COMPONENT lib FILES_MATCHING PATTERN "*.py")
install(FILES ${external_files} ${readme_file} ${CEPGEN_ADDONS_FILE}
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/CepGen
  COMPONENT lib)
export(PACKAGE CepGen)
include(BuildPackage)

#----- documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message(STATUS "Doxygen found")
  add_custom_target(doc_doxygen
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Doxygen" VERBATIM)
endif()
