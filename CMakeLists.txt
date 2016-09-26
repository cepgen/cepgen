cmake_minimum_required(VERSION 2.4)
project(CepGen)
set(PROJECT_VERSION 1)

#----- enable fortran (for Pythia/Herwig/...)
if(NOT CMAKE_GENERATOR STREQUAL Xcode)
  enable_language(Fortran OPTIONAL)
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -gfortran")

include_directories(external)
include_directories(export)
