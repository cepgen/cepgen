#--- searching for Photos++
find_library(PHOTOSPP Photospp HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_library(PHOTOSPP_HepMC3 PhotosppHepMC3 HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_path(PHOTOSPP_INCLUDE Photos HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES include)

if(NOT PHOTOSPP OR NOT PHOTOSPP_HepMC3)
  return()
endif()

#----- build the object
cepgen_build(CepGenPhotos
    SOURCES src/*.cpp
    LIBRARIES ${PHOTOSPP} ${PHOTOSPP_HepMC3} CepGenHepMC3
    INCLUDES ${PHOTOSPP_INCLUDE}
    COMPONENT photos)
cpack_add_component(photos
    DISPLAY_NAME "CepGen Photos++ wrapper library"
    DESCRIPTION "Collection of CepGen wrappers to the Photos library"
    DEPENDS lib)
