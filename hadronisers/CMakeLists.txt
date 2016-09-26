file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- enable fortran (for Pythia/Herwig/...)
if(NOT CMAKE_GENERATOR STREQUAL Xcode)
  enable_language(Fortran OPTIONAL)
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/external/pythia6*.f")
  message(STATUS "Pythia6 found!")
  add_definitions(-DPYTHIA6)
else()
  message(STATUS "Pythia6 not found. Will proceed without it")
endif()

add_library(Hadronisers SHARED ${sources})
