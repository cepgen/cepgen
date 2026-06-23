#--- searching for PGS
set(PGS_DIRS $ENV{PGS_DIR} ${PGS_DIR})
find_library(PGS pgslib HINTS ${PGS_DIRS} PATH_SUFFIXES lib)
find_library(PGS_STDHEP stdhep HINTS ${PGS_DIRS} PATH_SUFFIXES lib)
if(NOT PGS OR NOT PGS_STDHEP)
  return()
endif()

#----- build the object

cepgen_build(CepGenPGS SOURCES src/*.cpp
    EXT_LIBS ${PGS} ${PGS_STDHEP}
    INSTALL_COMPONENT pgs)
