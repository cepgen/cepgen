#--- searching for Tauola++
find_library(TAUOLAPP NAMES TauolaCxxInterface HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES lib)
find_library(TAUOLAPP_Fortran NAMES TauolaFortran HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES lib)
find_library(TAUOLAPP_HepMC3 TauolaHepMC3 HINTS $ENV{TAUOLAPP_DIR} PATH_SUFFIXES lib)
find_path(TAUOLAPP_INCLUDE Tauola HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES include)

#--- searching for HepMC3 headers
set(HEPMC3_DIRS $ENV{HEPMC3_DIR} /usr /usr/local)
find_path(HEPMC3_INCLUDE NAMES HepMC3 HepMC HINTS ${HEPMC3_DIRS} PATH_SUFFIXES include)

if(NOT TAUOLAPP OR NOT TAUOLAPP_HepMC3)  # Fortran standalone library is left optional
  return()
endif()

#----- build the object
cepgen_build(CepGenTauola
    SOURCES src/TauolaFilter.cpp
    TESTS test/*.cc
    LIBRARIES ${TAUOLAPP} ${TAUOLAPP_Fortran} ${TAUOLAPP_HepMC3} CepGenHepMC3
    INCLUDES ${TAUOLAPP_INCLUDE} ${HEPMC3_INCLUDE}
    COMPONENT tauola)
cpack_add_component(tauola
    DISPLAY_NAME "CepGen Tauola++ wrapper library"
    DESCRIPTION "Collection of CepGen wrappers to the Tauola library"
    DEPENDS lib)
