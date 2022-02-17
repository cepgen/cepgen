#--- searching for HepMC3
set(HEPMC3_DIRS $ENV{HEPMC3_DIR} /usr /usr/local)
find_library(HEPMC3_LIB NAME HepMC3 HINTS ${HEPMC3_DIRS} PATH_SUFFIXES lib64 lib)
find_library(HEPMC3_ROOT_LIB NAME HepMC3rootIO HINTS ${HEPMC3_DIRS} PATH_SUFFIXES root)
find_path(HEPMC3_INCLUDE NAME HepMC3 HINTS ${HEPMC3_DIRS} PATH_SUFFIXES include)
find_path(HEPMC3_EXTRA_PLUGINS NAMES ConvertExample HINTS ${HEPMC3_DIRS} $ENV{HEPMC3_EXTRA_PLUGINS_DIR})

if(NOT HEPMC3_LIB)
  return()
endif()

set(libs ${HEPMC3_LIB})
set(defs)
#message(STATUS "HepMC version 3 found")
if(HEPMC3_ROOT_LIB)
  #message(STATUS "HepMC ROOT I/O library found")
  list(APPEND libs ${HEPMC3_ROOT_LIB})
  list(APPEND defs -DHEPMC3_ROOTIO)
endif()

#----- build the object

cepgen_build(CepGenHepMC3 SOURCES *.cpp
    EXT_LIBS ${libs}
    EXT_HEADERS ${HEPMC3_INCLUDE}
    DEFINITIONS ${defs}
    INSTALL_COMPONENT hepmc3)
