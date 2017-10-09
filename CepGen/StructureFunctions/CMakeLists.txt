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

set(SCH_PATH ${PROJECT_SOURCE_DIR}/External/F2_Luxlike_fit)
file(GLOB sch_sources ${SCH_PATH}/F2_luxlike_minimal.f
                      ${SCH_PATH}/F2_perturbative.f ${SCH_PATH}/alphaS.f ${SCH_PATH}/mstw*.f
                      ${SCH_PATH}/cepgen_f2.cpp)
if(sch_sources)
  message(STATUS "LUXlike structure functions found!")
  add_definitions(-DSchaeferF2)
  list(APPEND sources ${sch_sources})
  configure_file(${SCH_PATH}/mstw_f2_scan_nnlo.txt ${PROJECT_BINARY_DIR}/External/F2_Luxlike_fit/mstw_f2_scan_nnlo.txt COPYONLY)
else()
  message(STATUS "LUXlike structure functions not found")
endif()

#----- build the object

add_library(CepGenStructureFunctions SHARED ${sources})
target_link_libraries(CepGenStructureFunctions ${external_deps})

#----- installation rules

install(TARGETS CepGenStructureFunctions DESTINATION lib)

