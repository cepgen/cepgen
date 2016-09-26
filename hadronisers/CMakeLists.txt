file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- enable fortran (for Pythia/Herwig/...)
if(NOT CMAKE_GENERATOR STREQUAL Xcode)
  enable_language(Fortran OPTIONAL)
endif()

#----- check the external dependencies

if(EXISTS "${PROJECT_SOURCE_DIR}/external/pythia6*.f")
  message(STATUS "Pythia6 found!")
  add_definitions(-DPYTHIA6)
else()
  message(STATUS "Pythia6 not found. Will proceed without it")
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/external/jetset7*.f")
  message(STATUS "Jetset7 found!")
  add_definitions(-DJETSET7)
else()
  message(STATUS "Jetset7 not found. Will proceed without it")
endif()

#----- build the object

add_library(Hadronisers OBJECT ${sources})
