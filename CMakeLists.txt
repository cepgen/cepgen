cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -g")

#----- define all individual modules to be built beforehand

set(MODULES core processes hadronisers physics export)

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

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
foreach(_module ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_module})
  add_subdirectory(${_module})
endforeach()

#----- link everything into the executables

add_executable(cepgen main.cpp)
target_link_libraries(cepgen CepGenCore CepGenExternalHadronisers CepGenProcesses)

find_package(HepMC)
if(HEPMC_LIBRARIES)
  add_executable(cepgen-lhe cepgen-lhe.cpp)
  target_link_libraries(cepgen-lhe CepGenCore CepGenExternalHadronisers CepGenProcesses CepGenExporter)
endif()

#----- copy the input cards and other files

file(GLOB input_cards RELATIVE ${PROJECT_SOURCE_DIR} cards/*)
foreach(_files ${input_cards})
  configure_file(
    ${_files}
    ${_files}
    COPYONLY
)
endforeach()

file(GLOB readme RELATIVE ${PROJECT_SOURCE_DIR} README)
configure_file(${readme} ${readme} COPYONLY)

#----- set the tests directory

add_subdirectory("test")
