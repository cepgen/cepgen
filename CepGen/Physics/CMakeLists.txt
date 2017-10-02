file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenPhysics SHARED ${sources})

#----- check the external dependencies

set(external_deps "")
if(EXISTS $ENV{has_grv})
  add_library(GRV_PDF SHARED $ENV{has_grv})
  list(APPEND ${external_deps} GRV_PDF)
  install(TARGETS GRV_PDF DESTINATION lib)
endif()
if(EXISTS $ENV{has_schaefer})
  add_library(SchaeferF2 SHARED ${sources})
  list(APPEND ${external_deps} SchaeferF2)
  install(TARGETS SchaeferF2 DESTINATION lib)
endif()
target_link_libraries(CepGenPhysics ${external_deps})

install(TARGETS CepGenPhysics DESTINATION lib)

