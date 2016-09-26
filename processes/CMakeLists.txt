file(GLOB sources *.cpp)
file(GLOB headers *.h)

include_directories(${PROJECT_SOURCE_DIR})

add_library(Processes SHARED ${sources})

file(GLOB has_grv ${PROJECT_SOURCE_DIR}/external/grv_*.f)
if(NOT ${has_grv})
  message(STATUS "GRV PDFset found in ${has_grv}!")
  add_definitions(-DGRVPDF)
  add_library(grv_shared STATIC ${has_grv})
  target_link_libraries(Processes grv_shared)
else()
  message(STATUS "GRV PDFset not found. Will proceed without it")
endif()
