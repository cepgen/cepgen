file(GLOB sf_sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenStructureFunctions SHARED ${sf_sources})

#----- check the external dependencies

set(external_deps "")
if(EXISTS $ENV{has_grv})
  add_library(GRV_PDF SHARED $ENV{has_grv})
  list(APPEND external_deps GRV_PDF)
  install(TARGETS GRV_PDF DESTINATION lib)
endif()
if(EXISTS $ENV{has_schaefer})
  list(APPEND external_deps SchaeferF2)
endif()
target_link_libraries(CepGenStructureFunctions ${external_deps})

#----- installation rules

install(TARGETS CepGenStructureFunctions DESTINATION lib)

