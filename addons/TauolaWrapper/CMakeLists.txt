#--- searching for Tauola++
find_library(TAUOLAPP NAMES TauolaCxxInterface HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES lib lib64)
find_library(TAUOLAPP_Fortran NAMES TauolaFortran HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES lib lib64)
find_library(TAUOLAPP_HepMC3 TauolaHepMC3 HINTS $ENV{TAUOLAPP_DIR} PATH_SUFFIXES lib lib64)
find_path(TAUOLAPP_INCLUDE Tauola HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES include)

if (NOT TAUOLAPP OR NOT TAUOLAPP_HepMC3)  # Fortran standalone library is left optional
    return()
endif ()

#----- build the object
cepgen_build(CepGenTauola
        SOURCES src/TauolaFilter.cpp
        TESTS test/*.cc
        LIBRARIES ${TAUOLAPP} ${TAUOLAPP_Fortran} ${TAUOLAPP_HepMC3} CepGenHepMC3
        INCLUDES ${TAUOLAPP_INCLUDE}
        COMPONENT tauola)
cpack_add_component(tauola
        DISPLAY_NAME "CepGen Tauola++ wrapper library"
        DESCRIPTION "Collection of CepGen wrappers to the Tauola library"
        DEPENDS lib)
