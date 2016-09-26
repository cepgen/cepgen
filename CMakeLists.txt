cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CepGen)
set(PROJECT_VERSION 1)

#----- define all individual modules to be built beforehand

set(MODULES core processes hadronisers physics export)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -gfortran")

#----- build all the intermediate objects

foreach(_module ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_module})
  add_subdirectory(${_module})
endforeach()

#----- link everything into the executables

add_executable(cepgen main.cpp $<TARGET_OBJECTS:Core> $<TARGET_OBJECTS:Physics> $<TARGET_OBJECTS:Processes> $<TARGET_OBJECTS:Hadronisers>)
