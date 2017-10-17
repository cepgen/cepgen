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
set(CEPGEN_LIBRARIES CepGenCore CepGenEvent CepGenPhysics CepGenProcesses CepGenStructureFunctions CepGenExporter CepGenExternalHadronisers)

#----- define all individual modules to be built beforehand

set(CEPGEN_MODULES Processes Event Physics StructureFunctions Export Hadronisers)

#----- enable fortran (for Pythia/Herwig/...)

enable_language(Fortran)

#----- include external paths

include(UseEnvironment)
set(CEPGEN_EXTERNAL_REQS "")
if(GSL_LIB)
  message(STATUS "GSL found in ${GSL_LIB}")
  list(APPEND CEPGEN_EXTERNAL_REQS ${GSL_LIB} ${GSL_CBLAS_LIB})
else()
  message(FATAL_ERROR "GSL not found on your system!")
endif()
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
  add_definitions(-DLIBLHAPDF)
  list(APPEND CEPGEN_EXTERNAL_REQS ${LHAPDF})
  include_directories(${LHAPDF_INCLUDE})
endif()
if(HEPMC_LIB)
  message(STATUS "HepMC found in ${HEPMC_INCLUDE}")
  include_directories(${HEPMC_INCLUDE})
  list(APPEND CEPGEN_EXTERNAL_REQS ${HEPMC_LIB})
  #list(APPEND CEPGEN_EXTERNAL_REQS ${HEPMC_FIO_LIB})
  add_definitions(-DLIBHEPMC)
endif()
if(PYTHIA8)
  message(STATUS "Pythia8 found in ${PYTHIA8}")
  list(APPEND CEPGEN_EXTERNAL_REQS ${PYTHIA8})
  include_directories(${PYTHIA8_INCLUDE})
endif()

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

