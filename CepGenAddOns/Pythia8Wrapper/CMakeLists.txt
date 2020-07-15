if(PYTHIA8)
  message(STATUS "Pythia 8 found in ${PYTHIA8}")
else()
  message(FATAL_ERROR "Trying to build Pythia8 wrapper library while no Pythia8 installation was found on system!")
endif()

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB sources *.cpp)

#----- build the objects

add_library(CepGenPythia8 SHARED ${sources})
target_link_libraries(CepGenPythia8 PRIVATE ${PYTHIA8} dl)
target_include_directories(CepGenPythia8 PRIVATE ${PYTHIA8_INCLUDE})
set_property(TARGET CepGenPythia8 PROPERTY CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-misleading-indentation")
list(APPEND CEPGEN_ADDONS_EXT CepGenPythia8)
set(CEPGEN_ADDONS_EXT ${CEPGEN_ADDONS_EXT} PARENT_SCOPE)

#----- installation rules

install(TARGETS CepGenPythia8 DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT pythia8)
install(DIRECTORY DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT pythia8 FILES_MATCHING PATTERN "*.h")