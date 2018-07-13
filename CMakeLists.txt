cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

#----- include external paths

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake)
set(CEPGEN_EXTERNAL_CORE_REQS "")
set(CEPGEN_EXTERNAL_IO_REQS "")
set(CEPGEN_EXTERNAL_HADR_REQS "")
set(CEPGEN_EXTERNAL_STRF_REQS "")
include(UseEnvironment)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_LIBRARIES CepGenCore CepGenEvent CepGenPhysics CepGenProcesses CepGenStructureFunctions CepGenIO CepGenExternalHadronisers)

#----- define all individual modules to be built beforehand

set(CEPGEN_MODULES Processes Event Physics StructureFunctions IO Hadronisers)

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

#----- installation rules

set(MODS "")
foreach(_module ${CEPGEN_MODULES})
  list(APPEND MODS CepGen/${module})
endforeach()
install(DIRECTORY ${MODS} DESTINATION include/CepGen FILES_MATCHING PATTERN "*.h")

#----- set the tests/utils directory

add_subdirectory(test)

