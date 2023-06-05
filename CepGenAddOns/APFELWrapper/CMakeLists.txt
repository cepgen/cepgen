#--- searching for APFEL
find_library(APFEL APFEL HINTS $ENV{APFEL_DIR} PATH_SUFFIXES lib)
find_path(APFEL_INCLUDE APFEL HINTS $ENV{APFEL_DIR} PATH_SUFFIXES include)

if(NOT APFEL)
  return()
endif()

#----- build the object

cepgen_build(CepGenAPFEL SOURCES AlphaSAPFEL.cpp PartonicStructureFunctionsAPFEL.cpp
    EXT_LIBS ${APFEL}
    EXT_HEADERS ${APFEL_INCLUDE}
    INSTALL_COMPONENT apfel)
cpack_add_component(apfel
    DISPLAY_NAME "CepGen APFEL wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the APFEL library"
    DEPENDS lib)
