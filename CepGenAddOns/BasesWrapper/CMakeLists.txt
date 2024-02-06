#--- searching for Bases
find_library(BASES bases HINTS $ENV{BASES_DIR})
if(NOT BASES)
  return()
endif()

#----- build the object
cepgen_build(CepGenBases SOURCES *.cpp
    EXT_LIBS ${BASES}
    PROPERTY POSITION_INDEPENDENT_CODE ON
    INSTALL_COMPONENT bases)
