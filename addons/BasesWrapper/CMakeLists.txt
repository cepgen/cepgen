#--- searching for Bases
find_library(BASES bases HINTS $ENV{BASES_DIR})
if(NOT BASES)
  return()
endif()

#----- build the object
cepgen_build(CepGenBases
    SOURCES src/*.cpp
    LIBRARIES ${BASES}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
    PROPERTY POSITION_INDEPENDENT_CODE ON
    COMPONENT bases)
install(DIRECTORY CepGenBases
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    COMPONENT bases)
