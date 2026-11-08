#--- searching for REvolver
find_library(REVOLVER REvolver HINTS $ENV{REVOLVER_DIR} PATH_SUFFIXES lib build)
find_path(REVOLVER_INCLUDE REvolver.h HINTS $ENV{REVOLVER_DIR} PATH_SUFFIXES include code/include)

if(NOT REVOLVER)
  return()
endif()

#----- build the object
cepgen_build(CepGenREvolver
    SOURCES src/*.cpp
    LIBRARIES ${REVOLVER}
    INCLUDES ${REVOLVER_INCLUDE}
    COMPONENT revolver)
