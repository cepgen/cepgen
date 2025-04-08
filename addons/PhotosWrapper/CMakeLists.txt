#--- searching for Photos++
find_library(PHOTOSPP Photospp HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_library(PHOTOSPP_HEPMC3 PhotosppHepMC3 HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_path(PHOTOSPP_INCLUDE Photos HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES include)

#--- searching for HepMC3 headers
set(HEPMC3_DIRS $ENV{HEPMC3_DIR} /usr /usr/local)
find_path(HEPMC3_INCLUDE NAMES HepMC3 HepMC HINTS ${HEPMC3_DIRS} PATH_SUFFIXES include)

if(NOT PHOTOSPP)
  return()
endif()

#----- build the object
cepgen_build(CepGenPhotos
    SOURCES src/*.cpp
    LIBRARIES ${PHOTOSPP} ${PHOTOSPP_HEPMC3} CepGenHepMC3
    INCLUDES ${PHOTOSPP_INCLUDE} ${HEPMC3_INCLUDE}
    COMPONENT photos)
cpack_add_component(photos
    DISPLAY_NAME "CepGen Photos++ wrapper library"
    DESCRIPTION "Collection of CepGen wrappers to the Photos library"
    DEPENDS lib)
