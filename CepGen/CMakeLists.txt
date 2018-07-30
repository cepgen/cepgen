file(GLOB core_sources Core/*.cpp *.cpp)
file(GLOB phys_sources Physics/*.cpp)
file(GLOB sf_sources StructureFunctions/*.cpp)

include_directories(${PROJECT_SOURCE_DIR})

#----- check the external dependencies for SFs

set(GRV_PATH ${PROJECT_SOURCE_DIR}/External)
file(GLOB grv_sources ${GRV_PATH}/grv_*.f)
if(grv_sources)
  message(STATUS "GRV PDFset found in ${grv_sources}!")
  add_definitions(-DGRVPDF)
  list(APPEND sf_sources ${grv_sources})
else()
  message(STATUS "GRV PDFset not found. Will proceed without it")
endif()

#----- check the external dependencies for physics utilities

if(alphas_sources)
  list(APPEND phys_sources ${alphas_sources})
endif()

#----- build the object

add_library(CepGenCore SHARED ${core_sources} ${phys_sources} ${sf_sources})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_CORE_REQS})
target_link_libraries(CepGenCore ${CEPGEN_EXTERNAL_STRF_REQS})

#----- installation rules

install(TARGETS CepGenCore DESTINATION lib)

