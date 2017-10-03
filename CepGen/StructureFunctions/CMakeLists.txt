file(GLOB sf_sources *.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- build the object

add_library(CepGenStructureFunctions SHARED ${sf_sources})

#----- check the external dependencies

if(EXISTS $ENV{has_grv})
  add_library(GRV_PDF SHARED $ENV{has_grv})
  target_link_libraries(CepGenStructureFunctions GRV_PDF)
endif()

install(TARGETS CepGenStructureFunctions DESTINATION lib)
install(TARGETS GRV_PDF DESTINATION lib)

