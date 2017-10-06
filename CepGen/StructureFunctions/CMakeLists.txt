file(GLOB sf_sources *.cpp)
set(sources ${sf_sources})

include_directories(${PROJECT_SOURCE_DIR})

#----- check the external dependencies

set(GRV_PATH ${PROJECT_SOURCE_DIR}/External)
file(GLOB grv_sources ${GRV_PATH}/grv_*.f)
if(grv_sources)
  message(STATUS "GRV PDFset found in ${grv_sources}!")
  add_definitions(-DGRVPDF)
  list(APPEND sources ${grv_sources})
else()
  message(STATUS "GRV PDFset not found. Will proceed without it")
endif()

#----- build the object

add_library(CepGenStructureFunctions SHARED ${sources})
target_link_libraries(CepGenStructureFunctions ${external_deps})

#----- installation rules

install(TARGETS CepGenStructureFunctions DESTINATION lib)

