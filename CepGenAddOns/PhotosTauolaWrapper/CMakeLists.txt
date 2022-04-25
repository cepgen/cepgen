#--- searching for Photos++
find_library(PHOTOSPP Photospp HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES lib)
find_path(PHOTOSPP_INCLUDE Photos HINTS $ENV{PHOTOSPP_DIR} PATH_SUFFIXES include)

if(NOT PHOTOSPP OR NOT PHOTOSPP_INCLUDE)
  return()
endif()

#----- build the object

cepgen_build(CepGenPhotosTauola SOURCES *.cpp
    EXT_LIBS ${PHOTOSPP}
    EXT_HEADERS ${PHOTOSPP_INCLUDE}
    INSTALL_COMPONENT photostauola)
