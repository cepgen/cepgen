file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- enable fortran (for Pythia/Herwig/...)
if(NOT CMAKE_GENERATOR STREQUAL Xcode)
  enable_language(Fortran OPTIONAL)
endif()

#----- build the object

add_library(Hadronisers SHARED ${sources})

#----- check the external dependencies

file(GLOB has_pythia6 ${PROJECT_SOURCE_DIR}/external/pythia6*.f)
if(NOT ${has_pythia})
  message(STATUS "Pythia6 found in ${has_pythia}!")
  add_definitions(-DPYTHIA6)
else()
  message(STATUS "Pythia6 not found. Will proceed without it")
endif()

file(GLOB has_herwig ${PROJECT_SOURCE_DIR}/external/herwig_*.f)
if(NOT ${has_herwig})
  message(STATUS "Herwig6 found in ${has_herwig}!")
  add_definitions(-DHERWIG6)
else()
  message(STATUS "Herwig6 not found. Will proceed without it")
endif()

file(GLOB has_jetset ${PROJECT_SOURCE_DIR}/external/jetset*.f)
if(NOT ${has_jetset})
  message(STATUS "Jetset7 found in ${has_jetset}!")
  add_definitions(-DJETSET7)
else()
  message(STATUS "Jetset7 not found. Will proceed without it")
endif()

#target_link_libraries(Hadronisers 
