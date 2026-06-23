# search for YODA
set(YODA_DIRS $ENV{YODA_DIR} ${YODA_DIR} /usr /usr/local)
find_library(YODA YODA HINTS ${YODA_DIRS} PATH_SUFFIXES lib lib64)
find_path(YODA_INCLUDE YODA HINTS ${YODA_DIRS} PATH_SUFFIXES include)
if(NOT YODA)
  return()
endif()

# searching for Fastjet
find_library(FASTJET NAMES fastjet HINTS $ENV{FASTJET_DIR} ${FASTJET_DIR} PATH_SUFFIXES lib lib64)
find_path(FASTJET_INCLUDE NAMES fastjet HINTS $ENV{FASTJET_DIR} ${FASTJET_DIR} PATH_SUFFIXES include)
find_library(FASTJET_CONTRIB NAMES fastjetcontribfragile HINTS $ENV{FASTJET_CONTRIB_DIR} ${FASTJET_CONTRIB_DIR} PATH_SUFFIXES lib lib64)
if(NOT FASTJET)
  return()
endif()

# searching for Rivet
find_library(RIVET NAMES Rivet HINTS $ENV{RIVET_DIR} ${RIVET_DIR} PATH_SUFFIXES lib lib64)
find_path(RIVET_INCLUDE NAMES Rivet HINTS $ENV{RIVET_DIR} ${RIVET_DIR} PATH_SUFFIXES include)
if(NOT RIVET)
  return()
endif()

# transform the "a.b.c" YODA version string into an integer equivalent
find_file(YODA_VERSION_FILE YODA/Config/YodaConfig.h HINTS ${YODA_INCLUDE})
file(READ ${YODA_VERSION_FILE} version_file)
string(REGEX MATCH "YODA_VERSION \"([0-9]*)\.([0-9]*)\.([0-9]*)\"" _ ${version_file})
math(EXPR YODA_VERSION "${CMAKE_MATCH_1} * 10000 + ${CMAKE_MATCH_2} * 100 + ${CMAKE_MATCH_3}")

# build the object
cepgen_build(CepGenRivet
    SOURCES src/*.cpp
    LIBRARIES ${YODA} ${FASTJET} ${FASTJET_CONTRIB} ${RIVET} CepGenHepMC3
    INCLUDES ${YODA_INCLUDE} ${FASTJET_INCLUDE} ${RIVET_INCLUDE}
    TESTS test/*.cc
    DEFINITIONS -DYODA_VERSION=${YODA_VERSION}
    OPTIONS "-Wno-deprecated-copy"
    COMPONENT rivet)
cpack_add_component(rivet
    DISPLAY_NAME "CepGen Rivet wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Rivet library"
    DEPENDS lib)
