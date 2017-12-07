cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic-errors -std=c++11 -g")
else()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.7)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic-errors -std=c++11 -g")
  else()
    message(FATAL_ERROR "gcc version >= 4.7 is required")
  endif()
endif()
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/cmake)
set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_LIBRARIES CepGenCore CepGenEvent CepGenPhysics CepGenProcesses CepGenStructureFunctions CepGenExporter CepGenIO CepGenExternalHadronisers)

#----- define all individual modules to be built beforehand

set(CEPGEN_MODULES Processes Event Physics StructureFunctions Export IO Hadronisers)

#----- enable fortran (for Pythia/Herwig/...)

enable_language(Fortran)

#----- include external paths

set(CEPGEN_EXTERNAL_REQS "")
include(UseEnvironment)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(${CEPGEN_SOURCE_DIR})
foreach(_module ${CEPGEN_MODULES})
  include_directories(${CEPGEN_SOURCE_DIR}/${_module})
  add_subdirectory(${CEPGEN_SOURCE_DIR}/${_module})
endforeach()

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

