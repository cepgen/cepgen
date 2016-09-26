cmake_minimum_required(VERSION 2.4)
project(CepGen)
set(PROJECT_VERSION 1)
set(MODULES core hadronisers physics export)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -gfortran")

foreach(_module ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_module}/include)
  add_subdirectory(${_module})
endforeach()

