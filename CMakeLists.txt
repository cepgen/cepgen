cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

#SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.7)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic-errors -std=c++11 -g")
else()
  message(FATAL_ERROR "gcc version >= 4.7 is required")
endif()
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake)
set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_LIBRARIES CepGenCore CepGenPhysics CepGenProcesses)

#----- define all individual modules to be built beforehand

set(CEPGEN_MODULES Processes Physics Export)

#----- enable fortran (for Pythia/Herwig/...)

enable_language(Fortran OPTIONAL)
add_subdirectory(External)

#----- include external paths

include(UseEnvironment)

set(CEPGEN_EXTERNAL_REQS ${GSL_LIB} ${GSL_CBLAS_LIB} ${CEPGEN_LIBRARIES})

if(MUPARSER)
  message(STATUS "muParser found in ${MUPARSER}")
  list(APPEND CEPGEN_EXTERNAL_REQS ${MUPARSER})
endif()
if(LIBCONFIG)
  message(STATUS "libconfig++ found in ${LIBCONFIG}")
  list(APPEND CEPGEN_EXTERNAL_REQS ${LIBCONFIG})
endif()
if(LHAPDF)
  message(STATUS "LHAPDF found in ${LHAPDF}")
  list(APPEND CEPGEN_EXTERNAL_REQS ${LHAPDF})
endif()
if(EXISTS $ENV{has_schaefer})
  message(STATUS "Schaefer's F2 found!")
  list(APPEND CEPGEN_EXTERNAL_REQS SchaeferF2)
endif()

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(${CEPGEN_SOURCE_DIR})
foreach(_module ${CEPGEN_MODULES})
  include_directories(${CEPGEN_SOURCE_DIR}/${_module})
  add_subdirectory(${CEPGEN_SOURCE_DIR}/${_module})
endforeach()

#----- optionally include the exports utils if HepMC is found

if(HEPMC_LIB)
  message(STATUS "HepMC found in ${HEPMC_INCLUDE}")
  include_directories(${HEPMC_INCLUDE})
  list(APPEND CEPGEN_LIBRARIES CepGenExporter)
endif()

#----- copy the input cards and other files

file(GLOB input_cards RELATIVE ${PROJECT_SOURCE_DIR} Cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${CEPGEN_SOURCE_DIR}/README README COPYONLY)

#----- installation rules

set(MODS "")
foreach(_module ${CEPGEN_MODULES})
  list(APPEND MODS CepGen/${module})
endforeach()
install(DIRECTORY ${MODS} DESTINATION include/CepGen FILES_MATCHING PATTERN "*.h")

#----- set the tests/utils directory

add_subdirectory(test)

