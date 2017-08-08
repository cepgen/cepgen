cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++11 -g")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake)
set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_LIBRARIES CepGenCore CepGenExternalHadronisers CepGenPhysics CepGenProcesses)

#----- define all individual modules to be built beforehand

set(CEPGEN_MODULES Core Processes Hadronisers Physics Export)

#----- enable fortran (for Pythia/Herwig/...)

enable_language(Fortran OPTIONAL)
add_subdirectory(external)

#----- include external paths

include(UseEnvironment)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
foreach(_module ${CEPGEN_MODULES})
  include_directories(${CEPGEN_SOURCE_DIR}/${_module})
  add_subdirectory(${CEPGEN_SOURCE_DIR}/${_module})
endforeach()

#----- optionally include the exports utils if HepMC is found

if(HEPMC_LIB)
  message(STATUS "HepMC found in " ${HEPMC_INCLUDE})
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

