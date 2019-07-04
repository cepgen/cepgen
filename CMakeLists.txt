cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake $ENV{ROOTSYS}/cmake)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-pg") # for gprof
#set(CMAKE_CXX_FLAGS_DEBUG "-fprofile-arcs -ftest-coverage") # for gcov
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

#----- include external paths

set(CEPGEN_EXTERNAL_CORE_REQS "")
set(CEPGEN_EXTERNAL_STRF_REQS "")
include(UseEnvironment)

set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_LIBRARIES CepGenCore CepGenEvent CepGenProcesses CepGenAddOns CepGenCards)

#----- define all individual modules to be built beforehand

set(CEPGEN_MODULES Processes Event Cards)

#----- enable fortran for external libraries linking

enable_language(Fortran)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(${CEPGEN_SOURCE_DIR})
foreach(_module ${CEPGEN_MODULES})
  add_subdirectory(${CEPGEN_SOURCE_DIR}/${_module})
endforeach()

#----- copy the input cards and other files

file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} Cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${CEPGEN_SOURCE_DIR}/README README COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/External/mstw_sf_scan_nnlo.dat ${PROJECT_BINARY_DIR}/External/mstw_sf_scan_nnlo.dat COPYONLY)

#----- installation rules

set(MODS "")
foreach(_module ${CEPGEN_MODULES})
  list(APPEND MODS CepGen/${module})
endforeach()
install(DIRECTORY ${MODS} DESTINATION include/CepGen FILES_MATCHING PATTERN "*.h")

#----- set the tests/utils directory

add_subdirectory(test)

#----- documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message(STATUS "Doxygen build started")
  add_custom_target(doc_doxygen COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating documentation with Doxygen" VERBATIM)
endif()

