#--- searching for CERNLIB
set(CERNLIB_DIRS $ENV{CERNLIB_DIR} ${CERNLIB_DIR})
find_library(CERNLIB_pawlib pawlib HINTS ${CERNLIB_DIRS} PATH_SUFFIXES lib)
find_library(CERNLIB_packlib packlib HINTS ${CERNLIB_DIRS} PATH_SUFFIXES lib)
find_path(CERNLIB_INCLUDE NAMES cfortran/cfortran.h cfortran/lapack.h PATHS ${CERNLIB_DIRS} PATH_SUFFIXES include)

if(NOT CERNLIB_pawlib OR NOT CERNLIB_packlib)
  return()
endif()

#----- build the object

cepgen_build(CepGenPAW SOURCES *.cpp
    EXT_LIBS ${CERNLIB_pawlib} ${CERNLIB_packlib} shift
    EXT_HEADERS ${CERNLIB_INCLUDE}
    DEFINITIONS -Df2cFortran
    CFLAGS "-Wall -W -Wno-unused-function -Wno-unused-label -Wno-unused-parameter -Wwrite-strings"
    INSTALL_COMPONENT paw)
