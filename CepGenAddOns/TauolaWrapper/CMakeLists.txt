#--- searching for Tauola++
find_library(TAUOLA NAMES TauolaCxxInterface HINTS $ENV{TAUOLA_DIR} ${TAUOLA_DIR} PATH_SUFFIXES lib)
find_path(TAUOLA_INCLUDE Tauola HINTS $ENV{TAUOLA_DIR} ${TAUOLA_DIR} PATH_SUFFIXES include)

if(NOT TAUOLA)
  return()
endif()

#----- build the object

cepgen_build(CepGenTauola SOURCES TauolaFilter.cpp PhotosTauolaInterface.cpp
    TESTS test/*.cc
    EXT_LIBS ${TAUOLA}
    EXT_HEADERS ${TAUOLA_INCLUDE}
    INSTALL_COMPONENT tauola)
