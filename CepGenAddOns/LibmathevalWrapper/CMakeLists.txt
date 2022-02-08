#--- searching for libmatheval
find_library(LIBMATHEVAL NAMES matheval HINTS $ENV{LIBMATHEVAL_DIR} PATH_SUFFIXES lib)
find_path(LIBMATHEVAL_INCLUDE NAMES matheval.h HINTS $ENV{LIBMATHEVAL_DIR} PATH_SUFFIXES include)

if(NOT LIBMATHEVAL)
  return()
endif()

file(GLOB sources *.cpp)

#----- build the object

cepgen_build(CepGenLibmatheval SOURCES ${sources}
    EXT_LIBS ${LIBMATHEVAL}
    EXT_HEADERS ${LIBMATHEVAL_INCLUDE}
    INSTALL_COMPONENT libmatheval)
