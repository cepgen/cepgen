cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

#----- define all individual modules to be built beforehand

set(MODULES core processes hadronisers physics export)

#----- enable fortran (for Pythia/Herwig/...)

enable_language(Fortran OPTIONAL)

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
foreach(_module ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_module})
  add_subdirectory(${_module})
endforeach()

#----- link everything into the executables

add_executable(cepgen main.cpp)
target_link_libraries(cepgen Core Physics Hadronisers Processes)

add_executable(cepgen-lhe cepgen-lhe.cpp)
target_link_libraries(cepgen-lhe Core Physics Hadronisers Processes Exporter)

#----- copy the input cards and other files

file(GLOB input_cards RELATIVE ${PROJECT_SOURCE_DIR} cards/*)
file(GLOB readme RELATIVE ${PROJECT_SOURCE_DIR} README)

foreach(_files ${input_cards})
  configure_file(
    ${_files}
    ${_files}
    COPYONLY
)
endforeach()
configure_file(${readme} ${readme} COPYONLY)
