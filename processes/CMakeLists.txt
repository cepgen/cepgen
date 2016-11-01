file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/physics)

add_library(CepGenProcesses SHARED ${sources})

file(GLOB has_grv ${PROJECT_SOURCE_DIR}/external/grv_*.f)
if(EXISTS ${has_grv})
  message(STATUS "GRV PDFset found in ${has_grv}!")
  add_definitions(-DGRVPDF)
  add_library(GRV_PDF SHARED ${has_grv})
  target_link_libraries(CepGenProcesses CepGenPhysics GRV_PDF)
else()
  message(STATUS "GRV PDFset not found. Will proceed without it")
endif()
