#--- searching for HepMC
set(HEPMC_DIRS $ENV{HEPMC_DIR} /usr /usr/local)
find_library(HEPMC_LIB NAMES HepMC3 HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES lib64 lib)
find_library(HEPMC_ROOT_LIB NAMES HepMC3rootIO HepMCrootIO HINTS ${HEPMC_DIRS} PATH_SUFFIXES root)
find_path(HEPMC_INCLUDE NAMES HepMC3 HepMC HINTS ${HEPMC_DIRS} PATH_SUFFIXES include)
find_path(HEPMC_EXTRA_PLUGINS NAMES ConvertExample HINTS ${HEPMC_DIRS} $HEPMC_EXTRA_PLUGINS_DIR $ENV{HEPMC_EXTRA_PLUGINS_DIR})

if(NOT HEPMC_LIB)
  return()
endif()

set(sources *.cpp)
set(libs ${HEPMC_LIB})
set(defs)
if(HEPMC_LIB MATCHES ".*HepMC3.?.so")
  #message(STATUS "HepMC version 3 found")
  list(APPEND defs -DHEPMC3)
  if(HEPMC_ROOT_LIB)
    #message(STATUS "HepMC ROOT I/O library found")
    list(APPEND libs ${HEPMC_ROOT_LIB})
    list(APPEND defs -DHEPMC3_ROOTIO)
  endif()
endif()

#----- build the object

cepgen_build(CepGenHepMC SOURCES *.cpp
    EXT_LIBS ${libs}
    EXT_HEADERS ${HEPMC_INCLUDE}
    DEFINITIONS ${defs}
    INSTALL_COMPONENT hepmc)
