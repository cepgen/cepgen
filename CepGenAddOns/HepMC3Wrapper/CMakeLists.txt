#--- searching for HepMC3
set(HEPMC3_DIRS $ENV{HEPMC3_DIR} /usr /usr/local)
find_library(HEPMC3_LIB NAME HepMC3 HINTS ${HEPMC3_DIRS} PATH_SUFFIXES lib64 lib)
find_library(HEPMC3_ROOT_LIB NAME HepMC3rootIO HINTS ${HEPMC3_DIRS} PATH_SUFFIXES root)
find_path(HEPMC3_INCLUDE NAMES HepMC3 HepMC HINTS ${HEPMC3_DIRS} PATH_SUFFIXES include)
find_path(HEPMC3_EXTRA_PLUGINS NAMES ConvertExample HINTS ${HEPMC3_DIRS} $ENV{HEPMC3_EXTRA_PLUGINS_DIR})

if(NOT HEPMC3_LIB OR NOT HEPMC3_INCLUDE)
  return()
endif()

set(libs ${HEPMC3_LIB})
set(headers ${HEPMC3_INCLUDE})
set(defs)
if(HEPMC3_ROOT_LIB)
  message(STATUS "HepMC ROOT I/O library found")
  list(APPEND libs ${HEPMC3_ROOT_LIB})
  list(APPEND defs -DHEPMC3_ROOTIO)
endif()

find_package(ZLIB)
if(ZLIB_FOUND)
  list(APPEND headers ${ZLIB_INCLUDE_DIR})
  list(APPEND libs ${ZLIB_LIBRARIES})
  list(APPEND defs HEPMC3_USE_COMPRESSION)
  list(APPEND defs HEPMC3_Z_SUPPORT=1)
endif()

find_package(LibLZMA)
if(LIBLZMA_FOUND)
  list(APPEND headers ${LIBLZMA_INCLUDE_DIRS})
  list(APPEND libs ${LIBLZMA_LIBRARIES})
  list(APPEND defs HEPMC3_USE_COMPRESSION)
  list(APPEND defs HEPMC3_LZMA_SUPPORT=1)
endif()

find_package(BZip2)
if(BZIP2_FOUND)
  list(APPEND headers ${BZIP2_INCLUDE_DIRS})
  list(APPEND libs ${BZIP2_LIBRARIES})
  list(APPEND defs HEPMC3_USE_COMPRESSION)
  list(APPEND defs HEPMC3_BZ2_SUPPORT=1)
endif()

#----- build the object

cepgen_build(CepGenHepMC3 SOURCES *.cpp
    EXT_LIBS ${libs}
    EXT_HEADERS ${headers}
    DEFINITIONS ${defs}
    INSTALL_COMPONENT hepmc3)
cpack_add_component(hepmc3
    DISPLAY_NAME "CepGen HepMC3 wrappers library"
    DESCRIPTION "Collection of CepGen wrappers to the HepMC3 library"
    DEPENDS lib)
