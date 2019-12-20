set(CPACK_RPM_COMPONENT_INSTALL ON)
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)
set(VERSION 0.9.7)

option(CMAKE_BUILD_CORE "Build Core library" ON)
option(CMAKE_BUILD_ADDONS "Build AddOns" ON)
option(CMAKE_BUILD_TESTS "Build tests" OFF)

#----- release build by default

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#----- check that at least one target is to be built

if(NOT CMAKE_BUILD_CORE AND NOT CMAKE_BUILD_ADDONS)
  message(FATAL_ERROR "Nothing to be built!")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake $ENV{ROOTSYS}/cmake)
set(CMAKE_CXX_FLAGS_DEBUG "-pg") # for gprof
#set(CMAKE_CXX_FLAGS_DEBUG "-fprofile-arcs -ftest-coverage") # for gcov
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(LIB_INSTALL_DIR lib)

#----- include external paths

include(GNUInstallDirs)
include(UseEnvironment)
include(BuildRPM)

set(CEPGEN_LIBRARIES)
set(CEPGEN_LIB_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

#----- define all individual modules/libraries to be built beforehand

set(CEPGEN_MODULES Event Cards)

#----- enable fortran for external libraries linking

enable_language(Fortran)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
if(CMAKE_BUILD_CORE)
  add_subdirectory(CepGen)
  foreach(_module ${CEPGEN_MODULES})
    add_subdirectory(${CEPGEN_LIB_DIR}/${_module})
  endforeach()
  list(APPEND CEPGEN_LIBRARIES CepGenCore CepGenProcesses CepGenEvent CepGenCards)
endif()
if(CMAKE_BUILD_ADDONS)
  add_subdirectory(CepGenAddOns)
  list(APPEND CEPGEN_LIBRARIES CepGenAddOns)
endif()
export(PACKAGE cepgen)

#----- copy the input cards and other files

file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} Cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
file(GLOB external_files ${PROJECT_SOURCE_DIR}/External/mstw_sf_scan_nnlo.dat ${PROJECT_SOURCE_DIR}/External/*.mcd)
file(COPY ${external_files} DESTINATION ${PROJECT_BINARY_DIR}/External)
file(GLOB readme_file ${PROJECT_SOURCE_DIR}/CepGen/README)
file(COPY ${readme_file} DESTINATION ${PROJECT_BINARY_DIR})

#----- installation rules

set(MODS)
foreach(_module ${CEPGEN_MODULES})
  list(APPEND MODS CepGen/${_module})
endforeach()
foreach(_module ${CEPGEN_LIBRARIES})
  message(STATUS ">> ${_module}")
endforeach()
install(DIRECTORY ${MODS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/CepGen COMPONENT devel FILES_MATCHING PATTERN "*.h")

#----- set the tests/utils directory

if(CMAKE_BUILD_CORE)
  #--- build the executables
  file(GLOB_RECURSE exec_sources RELATIVE ${CEPGEN_SOURCE_DIR} ${CEPGEN_SOURCE_DIR}/*.cpp)
  foreach(exec_src ${exec_sources})
    string(REPLACE ".cpp" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${CEPGEN_SOURCE_DIR}/${exec_src})
    target_link_libraries(${exec_bin} ${CEPGEN_LIBRARIES})
  endforeach()
  #--- build the tests
  if(CMAKE_BUILD_TESTS)
    enable_testing()
    include(CTest)
    add_subdirectory(test)
  endif()
endif()

#----- documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/doc/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message(STATUS "Doxygen found")
  add_custom_target(doc_doxygen COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating documentation with Doxygen" VERBATIM)
endif()

