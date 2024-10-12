#--- searching for TMDlib
set(TMDlib_DIRS $ENV{TMDlib_DIR} /usr /usr/local)
find_library(TMDlib TMDlib HINTS $TMDlib_DIRS} PATH_SUFFIXES lib lib64)
find_path(TMDlib_INCLUDE tmdlib/TMDlib.h HINTS ${TMDlib_DIRS} PATH_SUFFIXES include)

if(NOT TMDlib)
  return()
endif()

#----- build the object
cepgen_build(CepGenTMDlib
    SOURCES src/*.cpp
    LIBRARIES ${TMDlib}
    INCLUDES ${TMDlib_INCLUDE} ${TMDlib_INCLUDE}/tmdlib
    COMPONENT tmdlib)
