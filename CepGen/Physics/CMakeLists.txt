file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenPhysics SHARED ${sources})

#----- check the external dependencies

if(EXISTS $ENV{has_grv})
  add_library(GRV_PDF SHARED $ENV{has_grv})
  target_link_libraries(CepGenPhysics GRV_PDF)
endif()
