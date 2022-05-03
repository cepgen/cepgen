#--- searching for YODA
find_library(YODA YODA HINTS $ENV{YODA_DIR} PATH_SUFFIXES lib)
find_path(YODA_INCLUDE YODA HINTS $ENV{YODA_DIR} PATH_SUFFIXES include)

if(NOT YODA)
  return()
endif()

#----- build the object

cepgen_build(CepGenRivet SOURCES *.cpp
    EXT_LIBS ${YODA}
    EXT_HEADERS ${YODA_INCLUDE}
    OPTIONS "-Wno-deprecated-copy"
    INSTALL_COMPONENT rivet)
