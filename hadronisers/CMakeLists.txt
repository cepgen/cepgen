file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(Hadronisers SHARED ${sources})

#----- check the external dependencies

file(GLOB has_pythia6 ${PROJECT_SOURCE_DIR}/external/pythia6*.f)
if(NOT ${has_pythia6})
  message(STATUS "Pythia6 found in ${has_pythia6}!")
  add_definitions(-DPYTHIA6)
  add_library(pythia6_shared SHARED ${has_pythia6})
  target_link_libraries(Hadronisers pythia6_shared)
else()
  message(STATUS "Pythia6 not found. Will proceed without it")
endif()

file(GLOB has_herwig ${PROJECT_SOURCE_DIR}/external/herwig_*.f)
if(NOT ${has_herwig})
  message(STATUS "Herwig6 found in ${has_herwig}!")
  add_definitions(-DHERWIG6)
  add_library(herwig_shared SHARED ${has_herwig})
  target_link_libraries(Hadronisers herwig_shared)
else()
  message(STATUS "Herwig6 not found. Will proceed without it")
endif()

file(GLOB has_jetset ${PROJECT_SOURCE_DIR}/external/jetset*.f)
if(NOT ${has_jetset})
  message(STATUS "Jetset7 found in ${has_jetset}!")
  add_definitions(-DJETSET7)
  add_library(jetset_shared SHARED ${has_jetset})
  target_link_libraries(Hadronisers jetset_shared)
else()
  message(STATUS "Jetset7 not found. Will proceed without it")
endif()

#target_link_libraries(Hadronisers 
