# search for YODA
set(YODA_DIRS $ENV{YODA_DIR} /usr /usr/local)
find_library(YODA YODA HINTS ${YODA_DIRS} PATH_SUFFIXES lib)
if(NOT YODA)
  return()
endif()
find_path(YODA_INCLUDE YODA HINTS ${YODA_DIRS} PATH_SUFFIXES include)

# transform the "a.b.c" YODA version string into an integer equivalent
find_file(YODA_VERSION_FILE YODA/Config/YodaConfig.h HINTS ${YODA_INCLUDE})
file(READ ${YODA_VERSION_FILE} version_file)
string(REGEX MATCH "YODA_VERSION \"([0-9]*)\.([0-9]*)\.([0-9]*)\"" _ ${version_file})
math(EXPR YODA_VERSION "${CMAKE_MATCH_1} * 10000 + ${CMAKE_MATCH_2} * 100 + ${CMAKE_MATCH_3}")

# build the object
cepgen_build(CepGenRivet SOURCES *.cpp
    EXT_LIBS ${YODA}
    EXT_HEADERS ${YODA_INCLUDE}
    DEFINITIONS -DYODA_VERSION=${YODA_VERSION}
    OPTIONS "-Wno-deprecated-copy"
    INSTALL_COMPONENT rivet)
cpack_add_component(rivet
    DISPLAY_NAME "CepGen Rivet wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the Rivet library"
    DEPENDS lib)
