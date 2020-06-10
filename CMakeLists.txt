set(CPACK_RPM_COMPONENT_INSTALL ON)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)
set(VERSION 0.9.7)

option(CMAKE_BUILD_CORE "Build Core library" ON)
option(CMAKE_BUILD_PROCESSES "Build physics processes" ON)
option(CMAKE_BUILD_ADDONS "Build AddOns" ON)
option(CMAKE_BUILD_SCAN "Build cross-section scan utilitary" OFF)
option(CMAKE_BUILD_CARD_CONV "Build steering cards converter utilitary" OFF)
option(CMAKE_BUILD_TESTS "Build tests" OFF)
option(CMAKE_BUILD_FOAM "Build FOAM integrator" OFF)
#----- release build by default

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#----- check that at least one target is to be built

if(NOT CMAKE_BUILD_CORE AND NOT CMAKE_BUILD_PROCESSES AND NOT CMAKE_BUILD_ADDONS)
  message(FATAL_ERROR "Nothing to be built!")
endif()

#set(CMAKE_INSTALL_RPATH "$ORIGIN")
#set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_PREFIX "/usr")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake $ENV{ROOTSYS}/cmake)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

set(CMAKE_CXX_FLAGS_DEBUG "-pg") # for gprof
#set(CMAKE_CXX_FLAGS_DEBUG "-fprofile-arcs -ftest-coverage") # for gcov
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -O2")
set(CMAKE_C_FLAGS_RELEASE "-O2")

#----- include external paths

include(GNUInstallDirs)
include(UseEnvironment)

set(CEPGEN_LIBRARIES)
set(CEPGEN_LIB_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_CORE_EXT ${GSL_LIB} ${GSL_CBLAS_LIB} dl)
set(CEPGEN_ADDONS_EXT)
set(CEPGEN_HAS_FUNCTIONALS OFF)

#----- enable fortran for external libraries linking

enable_language(Fortran)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
if(CMAKE_BUILD_PROCESSES)
  add_subdirectory(CepGenProcesses)
endif()
if(CMAKE_BUILD_ADDONS)
  add_subdirectory(CepGenAddOns)
  message(STATUS "Found addons: ${CEPGEN_ADDONS_EXT}")
endif()
if(CMAKE_BUILD_CORE)
  add_subdirectory(CepGen)
  list(APPEND CEPGEN_LIBRARIES CepGenCore)
  add_subdirectory(CepGen/Cards)
  list(APPEND CEPGEN_LIBRARIES CepGenCards)
endif()

#----- copy the input cards and other files

file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} Cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
file(GLOB external_files ${PROJECT_SOURCE_DIR}/External/mstw_sf_scan_nnlo.dat ${PROJECT_SOURCE_DIR}/External/*.mcd)
file(COPY ${external_files} DESTINATION ${PROJECT_BINARY_DIR})
file(GLOB readme_file ${PROJECT_SOURCE_DIR}/CepGen/README)
file(COPY ${readme_file} DESTINATION ${PROJECT_BINARY_DIR})

#----- set the tests/utils directory

if(CMAKE_BUILD_CORE)
  #--- build the executables
  add_subdirectory(src)
  #--- build the tests
  if(CMAKE_BUILD_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(test)
  endif()
endif()

#----- installation rules

install(DIRECTORY CepGen CepGenAddOns CepGenProcesses
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT devel FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Cards
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/CepGen
  COMPONENT lib FILES_MATCHING PATTERN "*.py")
install(FILES ${external_files} ${readme_file}
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
