cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++14 -g")
set(CEPGEN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/CepGen)
set(CEPGEN_LIBRARIES CepGenCore CepGenExternalHadronisers CepGenPhysics CepGenProcesses)

#----- define all individual modules to be built beforehand

set(MODULES Core Processes Hadronisers Physics Export)

#----- enable fortran (for Pythia/Herwig/...)

enable_language(Fortran OPTIONAL)

#----- include external hadronisers/...

add_subdirectory(external)

if(EXISTS $ENV{has_pythia6})
  add_definitions(-DPYTHIA6)
endif()
if(EXISTS $ENV{has_herwig})
  add_definitions(-DHERWIG)
endif()
if(EXISTS $ENV{has_jetset})
  add_definitions(-DJETSET)
endif()
if(EXISTS $ENV{has_grv})
  add_definitions(-DGRVPDF)
endif()

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
foreach(_module ${MODULES})
  include_directories(${CEPGEN_SOURCE_DIR}/${_module})
  add_subdirectory(${CEPGEN_SOURCE_DIR}/${_module})
endforeach()

#----- link everything into the executables

if($ENV{HOSTNAME} MATCHES "^lxplus[0-9]+.cern.ch")
  set(BASE_DIR "/cvmfs/sft.cern.ch/lcg/external")
  set(GSL_DIR "${BASE_DIR}/GSL/1.14/x86_64-slc5-gcc44-opt")
  set(HEPMC_DIR "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt")

  message(STATUS "Compiling on LXPLUS. Do not forget to source the environment variables!")

  find_library(GSL_LIB gsl HINTS "${GSL_DIR}/lib")
  find_library(GSL_CBLAS_LIB gslcblas HINTS "${GSL_DIR}/lib")
  #--- searching for HepMC
  find_library(HEPMC_LIB HepMC HINTS "${HEPMC_DIR}/lib")
  find_library(HEPMC_FIO_LIB HepMCfio HINTS "${HEPMC_DIR}/lib")
  find_path(HEPMC_INCLUDE HepMC HINTS "${HEPMC_DIR}/include")
else()
  find_library(GSL_LIB gsl)
  find_library(GSL_CBLAS_LIB gslcblas)
  find_library(HEPMC_LIB HepMC)
  find_path(HEPMC_INCLUDE HepMC)
endif()

if(HEPMC_LIB)
  message(STATUS "HepMC found in " ${HEPMC_INCLUDE})
  include_directories(${HEPMC_INCLUDE})
  list(APPEND CEPGEN_LIBRARIES CepGenExporter)
endif()

#----- copy the input cards and other files

file(GLOB input_cards RELATIVE ${PROJECT_SOURCE_DIR} cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()
configure_file(${CEPGEN_SOURCE_DIR}/README README COPYONLY)

#----- set the tests directory

add_subdirectory(test)
add_subdirectory(utils)
