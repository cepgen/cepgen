if(APFEL)
  message(STATUS "APFEL found in ${APFEL}")
else()
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources AlphaSAPFEL.cpp)

#----- build the objects

add_library(CepGenAPFEL SHARED ${sources})
target_link_libraries(CepGenAPFEL PRIVATE ${APFEL})
target_include_directories(CepGenAPFEL PRIVATE ${APFEL_INCLUDE})

#----- installation rules

install(TARGETS CepGenAPFEL DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT apfel)
