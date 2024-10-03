#--- searching for Pythia 8
set(PYTHIA8_DIRS $ENV{PYTHIA8_DIR} /usr /usr/local /opt/pythia8)
find_library(PYTHIA8 pythia8 HINTS ${PYTHIA8_DIRS} PATH_SUFFIXES lib lib64)
find_path(PYTHIA8_INCLUDE Pythia8 HINTS ${PYTHIA8_DIRS} PATH_SUFFIXES include/Pythia8 include/pythia8 include)

if(NOT PYTHIA8)
  return()
endif()

#--- searching for gzip
set(PYTHIA_DEFS)
find_program(GZIP_BIN gzip DOC "Path to gzip program")
if(GZIP_BIN)
  list(APPEND PYTHIA_DEFS "-DGZIP_BIN=\"${GZIP_BIN}\"")
endif()

#----- build the object
cepgen_build(CepGenPythia8 SOURCES src/*.cpp
    EXT_LIBS ${PYTHIA8} dl
    EXT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHIA8_INCLUDE}
    DEFINITIONS ${PYTHIA_DEFS}
    TESTS test/*.cc
    OPTIONS "-Wno-misleading-indentation"
    INSTALL_DIRECTORY CepGenPythia8
    INSTALL_COMPONENT pythia8)
cpack_add_component(pythia8
    DISPLAY_NAME "CepGen Pythia 8 wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to Pythia 8"
    DEPENDS lib)
