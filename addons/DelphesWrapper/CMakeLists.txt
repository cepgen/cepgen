#--- searching for ROOT
find_package(ROOT)
if(NOT ROOT_FOUND)
  return()
endif()

#--- searching for Delphes
find_library(DELPHES Delphes HINTS $ENV{DELPHES_DIR} PATH_SUFFIXES lib build)
if(NOT DELPHES)
  return()
endif()
find_path(DELPHES_INCLUDE NAMES modules classes HINTS $ENV{DELPHES_DIR} PATH_SUFFIXES include)
find_path(DELPHES_EXT_INCLUDE NAMES ExRootAnalysis HINTS $ENV{DELPHES_DIR} PATH_SUFFIXES external include)

#----- build the object
cepgen_build(CepGenDelphes SOURCES src/*.cpp
    EXT_HEADERS ${DELPHES_INCLUDE} ${DELPHES_EXT_INCLUDE} ${ROOT_INCLUDE_DIRS}
    EXT_LIBS ${DELPHES} ${ROOT_LIBRARIES}
    DEFINITIONS ${ROOT_DEFINITIONS}
    PROPERTY CXX_STANDARD ${ROOT_CXX_STANDARD}
    CFLAGS ${ROOT_CXX_FLAGS}
    INSTALL_COMPONENT delphes)
cpack_add_component(delphes
    DISPLAY_NAME "CepGen Delphes wrappers library"
    DESCRIPTION "A CepGen wrapper to the Delphes utility"
    DEPENDS lib)
