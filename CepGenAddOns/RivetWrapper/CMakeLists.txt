if(NOT YODA)
  message(FATAL_ERROR "Trying to build Rivet wrapper library while no Rivet installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources YODAHistsHandler.cpp)

#----- build the objects

add_library(CepGenRivet SHARED ${sources})
target_link_libraries(CepGenRivet PUBLIC ${YODA})
target_include_directories(CepGenRivet PUBLIC ${YODA_INCLUDE})

#----- installation rules

install(TARGETS CepGenRivet DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT rivet)
