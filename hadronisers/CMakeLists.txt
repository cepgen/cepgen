file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(Hadronisers SHARED ${sources})

#----- check the external dependencies

if(EXISTS $ENV{has_pythia6})
  add_library(pythia6_shared SHARED $ENV{has_pythia6})
  target_link_libraries(Hadronisers pythia6_shared)
endif()

if(EXISTS $ENV{has_herwig})
  add_library(herwig_shared SHARED $ENV{has_herwig})
  target_link_libraries(Hadronisers herwig_shared)
endif()

if(EXISTS $ENV{has_jetset})
  add_library(jetset_shared SHARED $ENV{has_jetset})
  target_link_libraries(Hadronisers jetset_shared)
endif()
