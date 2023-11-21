#--- searching for YODA
find_library(YODA YODA HINTS $ENV{YODA_DIR} ${YODA_DIR} PATH_SUFFIXES lib)
find_path(YODA_INCLUDE YODA HINTS $ENV{YODA_DIR} ${YODA_DIR} PATH_SUFFIXES include)
#--- searching for Fastjet
find_library(FASTJET NAMES fastjet HINTS $ENV{FASTJET_DIR} ${FASTJET_DIR} PATH_SUFFIXES lib)
find_path(FASTJET_INCLUDE NAMES fastjet HINTS $ENV{FASTJET_DIR} ${FASTJET_DIR} PATH_SUFFIXES include)
find_library(FASTJET_CONTRIB NAMES fastjetcontribfragile HINTS $ENV{FASTJET_CONTRIB_DIR} ${FASTJET_CONTRIB_DIR} PATH_SUFFIXES lib)
#--- searching for Rivet
find_library(RIVET NAMES Rivet HINTS $ENV{RIVET_DIR} ${RIVET_DIR} PATH_SUFFIXES lib)
find_path(RIVET_INCLUDE NAMES Rivet HINTS $ENV{RIVET_DIR} ${RIVET_DIR} PATH_SUFFIXES include)

if(NOT YODA OR NOT FASTJET OR NOT RIVET)
  return()
endif()

#----- build the object

cepgen_build(CepGenRivet SOURCES *.cpp
    EXT_LIBS ${RIVET} ${FASTJET} ${FASTJET_CONTRIB} ${YODA} CepGenHepMC2
    EXT_HEADERS ${RIVET_INCLUDE} ${FASTJET_INCLUDE} ${YODA_INCLUDE}
    TESTS test/*.cc
    OPTIONS "-Wno-deprecated-copy"
    INSTALL_COMPONENT rivet)
cpack_add_component(rivet
    DISPLAY_NAME "CepGen Rivet wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Rivet library"
    DEPENDS lib)
