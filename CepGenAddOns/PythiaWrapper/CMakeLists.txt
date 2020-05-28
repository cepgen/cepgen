if(PYTHIA8)
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
else()
  message(FATAL_ERROR "Trying to build Pythia8 wrapper library while no Pythia8 installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenPythia SHARED ${sources})
target_link_libraries(CepGenPythia PUBLIC ${PYTHIA8} dl)
target_include_directories(CepGenPythia PUBLIC ${PYTHIA8_INCLUDE})
list(APPEND CEPGEN_ADDONS_EXT CepGenPythia)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenPythia DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT pythia8)
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT pythia8 FILES_MATCHING PATTERN "*.h")
