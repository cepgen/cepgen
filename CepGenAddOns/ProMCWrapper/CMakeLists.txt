#--- linking with ProMC

set(PROMC_DIRS $ENV{PROMC_DIR} $ENV{PROMC} /usr /usr/local /usr/local/promc)
find_library(PROMC_LIB NAMES promc HINTS ${PROMC_DIRS} PATH_SUFFIXES lib)
find_path(PROMC_INCLUDE NAMES ProMCBook.h HINTS ${PROMC_DIRS} PATH_SUFFIXES src)
find_path(PROMC_EXT_INCLUDE NAMES CBook HINTS ${PROMC_DIRS} PATH_SUFFIXES include)

if(NOT PROMC_LIB)
  return()
endif()

#----- build the object

cepgen_build(CepGenProMC SOURCES *.cpp
    EXT_LIBS ${PROMC_LIB}
    EXT_HEADERS ${PROMC_INCLUDE} ${PROMC_EXT_INCLUDE}
    PROPERTY CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations"
    INSTALL_COMPONENT promc)
