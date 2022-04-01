#--- searching for LHAPDF
find_library(LHAPDF LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES lib)
find_path(LHAPDF_INCLUDE LHAPDF HINTS $ENV{LHAPDF_PATH} PATH_SUFFIXES include)

if(NOT LHAPDF)
  return()
endif()

#----- build the object

cepgen_build(CepGenLHAPDF SOURCES PartonicStructureFunctionsLHAPDF.cpp AlphaSLHAPDF.cpp
    EXT_LIBS ${LHAPDF}
    TESTS test/*.cc
    EXT_HEADERS ${LHAPDF_INCLUDE}
    INSTALL_COMPONENT lhapdf)
