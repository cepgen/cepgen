if(NOT PYTHIA8)
  message(FATAL_ERROR "Trying to build Pythia8 wrapper library while no Pythia8 installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenPythia SHARED ${sources})
target_link_libraries(CepGenPythia PUBLIC ${PYTHIA8} dl)
target_include_directories(CepGenPythia PUBLIC ${PYTHIA8_INCLUDE})

#----- installation rules

install(TARGETS CepGenPythia DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT pythia8)
