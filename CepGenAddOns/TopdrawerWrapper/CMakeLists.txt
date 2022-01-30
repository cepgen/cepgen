#--- searching for Topdrawer
find_program(TOPDRAWER_BIN NAMES td HINTS $ENV{TOPDRAWER_DIR} PATH_SUFFIXES bin)
#https://cp3.irmp.ucl.ac.be/projects/madgraph/raw-attachment/wiki/TopDrawer/td

if(NOT TOPDRAWER_BIN)
  return()
endif()

file(GLOB sources DrawerTopdrawer.cpp)

#----- build the object

cepgen_build(CepGenTopdrawer SOURCES ${sources}
    DEFINITIONS "-DTD_BIN=${TOPDRAWER_BIN}"
    INSTALL_COMPONENT topdrawer)
