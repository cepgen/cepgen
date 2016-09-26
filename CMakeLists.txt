cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

#----- define all individual modules to be built beforehand

set(MODULES core processes hadronisers physics export)

enable_language(Fortran)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

#----- build all the intermediate objects

include_directories(${PROJECT_SOURCE_DIR})
foreach(_module ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_module})
  add_subdirectory(${_module})
endforeach()

#----- link everything into the executables

#add_executable(cepgen main.cpp $<TARGET_OBJECTS:Core> $<TARGET_OBJECTS:Physics> $<TARGET_OBJECTS:Processes> $<TARGET_OBJECTS:Hadronisers>)
#target_link_libraries(cepgen ${GSL_LIBRARIES})

find_package(HepMC REQUIRED)

add_executable(cepgen-lhe cepgen-lhe.cpp $<TARGET_OBJECTS:Exporter>)
target_link_libraries(cepgen-lhe ${GSL_LIBRARIES} ${HEPMC_LIBRARIES} Core Physics Hadronisers Processes)

