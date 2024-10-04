#--- searching for Photos++
find_library(PHOTOSPP Photospp HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_library(PHOTOSPP_HEPMC3 PhotosppHepMC3 HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_path(PHOTOSPP_INCLUDE Photos HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES include)

#--- searching for Tauola++
find_library(TAUOLAPP NAMES TauolaCxxInterface HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES lib)
find_library(TAUOLAPP_HEPMC3 TauolaHepMC3 HINTS $ENV{TAUOLAPP_DIR} PATH_SUFFIXES lib)
find_path(TAUOLAPP_INCLUDE Tauola HINTS $ENV{TAUOLAPP_DIR} ${TAUOLAPP_DIR} PATH_SUFFIXES include)

#--- searching for HepMC3 headers
set(HEPMC3_DIRS $ENV{HEPMC3_DIR} /usr /usr/local)
find_path(HEPMC3_INCLUDE NAMES HepMC3 HepMC HINTS ${HEPMC3_DIRS} PATH_SUFFIXES include)

if(NOT PHOTOSPP AND NOT TAUOLAPP)
  return()
endif()

set(SRCS)
set(TESTS)
set(EXT_LIBS)
set(EXT_INCS)
if(PHOTOSPP)
  list(APPEND SRCS src/PhotosFilter.cpp)
  list(APPEND EXT_LIBS ${PHOTOSPP} ${PHOTOSPP_HEPMC3})
  list(APPEND EXT_INCS ${PHOTOSPP_INCLUDE})
endif()
if(TAUOLAPP)
  list(APPEND SRCS src/TauolaFilter.cpp)
  list(APPEND TESTS test/simple_tau_decay.cc)
  list(APPEND EXT_LIBS ${TAUOLAPP} ${TAUOLAPP_HEPMC3})
  list(APPEND EXT_INCS ${TAUOLAPP_INCLUDE})
endif()

#----- build the object
cepgen_build(CepGenPhotosTauola SOURCES ${SRCS}
    TESTS ${TESTS}
    EXT_LIBS ${EXT_LIBS} CepGenHepMC3
    EXT_HEADERS ${EXT_INCS} ${HEPMC3_INCLUDE}
    INSTALL_COMPONENT photostauola)
cpack_add_component(photostauola
    DISPLAY_NAME "CepGen Photos++/Tauola++ wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Photos and Tauola libraries"
    DEPENDS lib)
